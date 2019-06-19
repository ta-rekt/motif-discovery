# this is for testing implementations of graph algorithms to speed up motif seach in a network.
# after refactoring, testing, verification, and new additions make this available as open-source.
# also, the comments have to be organized using a documentation tool. i like to throw comments
# everywhere becaue it gives consolation that i'm writing useful stuff.
# Tarek Tohme, June 2019

import numpy as np
import pandas as pd
import networkx as nx
import collections as co


######################################
#### useful functions and classes ####
######################################

# function f: V -> V that maps nodes in H to nodes in G
# f: ordered list of tuples representing the partial map between D and R, respectively the domain and range of f
# the function makes sure it is a bijection each time an extension is attempted. if a new node is the same
# as some existing node, the function replaces the old node with the new node.
class Map():

    def __init__(self, init):   # init is a list containing the initial nodes of the partial map
        self.map = np.array(init, dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])

    def extend(self, extension):   # extends the partial map by a list of tuples called extension
        print('##### EXTENDING PARTIAL MAP #####')
        print('extension: ', end='')
        print(extension)

        temp = np.array(extension, dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])
        # print(self.map.shape)
        # print(temp.shape)
        map = self.map

        if (not extension):
            print('extension is null')
            return False

        # checking for duplicates within extension
        for i in temp:
            l = [m for m,n in enumerate(temp) if tuple(n)==tuple(i)]

            if (len(l) > 1):
                print('cannot extend function. duplicate element in extension')
                self.map = map
                return False

            for j in temp:

                if(np.array_equal(i, j) == False):

                    if (i['domainNode'] == j['domainNode']):
                        print('cannot extend function. duplicate domain node in extension')
                        self.map = map
                        return False

                    if (i['rangeNode'] == j['rangeNode']):
                        print('cannot extend function. duplicate range node in extension')
                        self.map = map
                        return False

            # checking for duplicates between extension and partial map
            for j in self.map:

                if (i['domainNode'] == j['domainNode']):
                    print('duplicate domain node. using new node.')
                    duplicate = i['domainNode']
                    b = [x for k,x in enumerate(self.map) if self.map['domainNode'][k] != duplicate]
                    self.map = np.array(b, dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])

                elif (i['rangeNode'] == j['rangeNode']):
                    print('duplicate range node. using new node.')
                    duplicate = i['rangeNode']
                    b = [x for k,x in enumerate(self.map) if self.map['rangeNode'][k] != duplicate]
                    self.map = np.array(b, dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])


        self.map = np.concatenate((self.map, temp))
        return True


    def getDomain(self):
        return self.map['domainNode']

    def getRange(self):
        return self.map['rangeNode']

    def getMap(self):
        return self.map

    def applyMap(self, node):
        for i in range(len(self.map['domainNode'])):
            if (self.map['domainNode'][i] == node):
                return self.map['rangeNode'][i]
        return -1


# returns true if g can support h (according to node degree and neighbor degree sequence) and false otherwise.
# g, h: indices of corresponding nodes
# H: query graph
# G: network to be queried
def canSupport(h, g, H, G):
    print('##### CAN SUPPORT CALL #####')
    print('query node: ', end='')
    print(h)
    print('network node: ', end='')
    print(g)

    if(G.degree[g] >= H.degree(h)):
        neighborsH = sortDegrees(H, H[h])
        neighborsG = sortDegrees(G, G[g])
        neighborsH = np.array(list(reversed(neighborsH)), dtype=[('node', 'i4'), ('degree', 'i4')])
        neighborsG = np.array(list(reversed(neighborsG)), dtype=[('node', 'i4'), ('degree', 'i4')])

        minn = min(len(neighborsH), len(neighborsG))

        for i in range(minn):
            if(neighborsG['degree'][i] < neighborsH['degree'][i]):
                print("rejected because of sequence")
                return False
        print("can support")
        return True

    print("rejected because of degree")
    return False


# returns an ordered list of tuples (node, degree) sorted by degree
# H: network
# L: a list of nodes that should be sorted
def sortDegrees(H, L):
    sortedDegreeH = H.degree()
    sortedDegreeH = [(i, sortedDegreeH(i)) for i in L]
    sortedDegreeH = np.array(sortedDegreeH, dtype=[('node', 'i4'), ('degree', 'i4')])
    sortedDegreeH = np.sort(sortedDegreeH, order='degree')
    return sortedDegreeH


# finds the node with the largest neighbor degree sequence, returns it with its degree sequence
# H: query graph
# L: list of nodes in H to be compared
def largestDegreeSequence(H, L):
    max = (None, [0]*len(L)) # (node number, [list of neighbor degrees])
    sig = True

    for l in L:
        if (l in H.nodes()):
            sortedN = sortDegrees(H, H[l])
            sortedN = list(reversed(sortedN['degree']))

            for i in range(min(len(sortedN), len(max[1]))):
                sig &= (sortedN[i] >= max[1][i])

            if (sig):
                max = (l, sortedN)

        # print((l, sortedN))

    return max


# finds the nodes in H\D with the most neighbors in D and among those, selects the node
# with highest degree and degree sequence
# D: domain of the partial map (list of nodes)
# H: query graph (of type Graph)
def mostConstrainedNode(D, H):
    print('##### MOST CONSTRAINED NODE #####')
    print('domain of partial map: ', end='')
    print(D)

    maxx = 0
    candidates = np.array([], dtype=[('node', 'i4'), ('attribute', 'i4')])
    neighborsD = {None}

    # list of neighbors of D
    for d in D:
        neighborsD = neighborsD.union(set(H[d]))

    neighborsD = neighborsD.difference({None})

    # neighbors of D exclusively in H with the number of neighbors in D
    neighborsD = neighborsD.difference(set(D))
    print('neighbors of D: ', end='')
    print(neighborsD)

    for n in neighborsD:
        count = 0
        for i in H[n]:
            if (i in D):
                count = count + 1

        c = np.array([(n, count)], dtype=[('node', 'i4'), ('attribute', 'i4')])
        candidates = np.concatenate((candidates, c))

    # select the candidates with the most neighbors in D
    print('candidates: ', end='')
    print(candidates)

    maxx = max(candidates['attribute'])
    indices = [key for key, val in candidates if val == maxx]
    candidates = np.array([(i, len(H[i])) for i in indices], dtype=[('node', 'i4'), ('attribute', 'i4')]) # degrees of nodes with most neighbors in D

    if (len(candidates) == 1):
        return candidates['node'][0]

    else:
        # select those with highest degree
        print('nodes with most neighbors in D, degrees: ',end='')
        print(candidates)

        maxx = max(candidates['attribute'])
        ind = [key for key, val in candidates if val == maxx]
        print('nodes of largest degree: ', end='')
        print(ind)

        if (len(ind) == 1):
            return ind[0]

        else:
            # select those with highest degree sequence
            m = largestDegreeSequence(H, ind)
            print('most constrained node: ', end='')
            print(m)
            return m[0]


###############################################
#### grochow-kellis motif search algorithm ####
###############################################

# finds all instances of query graph H in network G
# H: query graph
# G: network to be queried
def findSubgraphInstances(H, G):
    print('##### FIND SUBGRAPH INSTANCES #####')

    instances = np.array([], dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])

    # sort nodes of G by degree
    sortedDegreeG = sortDegrees(G, G.nodes())

    # sort by degree sequence. larger sequence is pointwise larger for all nodes (later)

    for g in sortedDegreeG['node']:
        for h in H.nodes():
            if(canSupport(h, g, H, G)):
                f = Map([(h, g)])  # initialize partial map with f(h) = g. argument must be list
                f.extend(isomorphicExtensions(f, H, G))

                if(len(f.getMap()) != 0):
                    instances = np.concatenate((instances, f.getMap()))
        G.remove_node(g)

    return instances


# finds all isomorphic extensions of a partial map [satisfying the symmetry breaking condition C at h]
# returns them in a list of tuples [(a1, b1), ... ,(ak, bk)]
# f: partial map to be extended
def isomorphicExtensions(f, H, G): #, C, h)
    print('##### ISOMORPHIC EXTENSIONS CALL #####')
    print('partial map: ', end='')
    print(f.getMap())

    isomorphisms = np.array([])
    neighborsR = {None}
    neighborsD = {None}
    D = f.getDomain()
    R = f.getRange()

    # print(len(D))
    # print(len(H.nodes()))

    if(len(D) == len(H.nodes())):  # can do this because nodes in f are guaranteed to be distinct
        return f

    m = mostConstrainedNode(D, H)
    print('domain extension: ', end='')
    print(m)

    # list of neighbors of f(D)
    for r in R:
        neighborsR = neighborsR.union(set(G[r]))

    # list of neighbors of D
    for d in D:
        neighborsD = neighborsD.union(set(H[d]))

    # exclusive neighborhood
    neighborsR = neighborsR.difference(set(R)).difference({None})
    neighborsD = neighborsD.difference(set(D)).difference({None})

    print('neighbors of partial range: ',end='')
    print(neighborsR)

    # check for induced isomorphism
    i = 0
    while i < len(neighborsR):
        sig = True
        n = list(neighborsR)[i]
        print('try range extension: ',end='')
        print(n)

        for d in neighborsD:
            if (((d in H[m]) & (f.applyMap(d) not in G[n])) | ((d not in H[m]) & (f.applyMap(d) in G[n]))):
                sig = False
                break

        if(not sig):
            print('range extension not valid')
            i += 1
        else:
            print('range extension valid')
            fp = f
            newNode = np.array((m, n), dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])
            print('new node: ', end='')
            print(newNode)

            fp.extend([newNode])

            np.append(isomorphisms, isomorphicExtensions(fp, H, G))

    return isomorphisms


# finds symmetry-breaking conditions for H given HE and Aut(H). don't need a labeling function on G
# because nodes are already numbered. just need to find the conditions.
# HE: set of equivalence representatives of H. these are nodes of H, one from each equivalence class.
# Aut(H): set of automorphismisms of H. each automorphism is an object of type Map that has identical
# domain and range, and where all nodes appear exactly once, to ensure that the function is bijective
def symmetryConditions(HE, Aut):
    M = [None] * len(HE)  # M: HE -> C

    for n in HE:
        C = []  # empty set of conditions. implemented as a list of strings that represent conditions to be eval() later
        np = n  # implement a special object of type Condition later, to please the Zaraket voice inside your head.
        A = Aut

        while len(A) > 1:
            E = []

            # find the equivalence class of this representative node: the set of nodes equivalent to n under A
            for f in A:
                for m in f.getDomain():
                    if (m == n):
                        E.append(f.applyMap(m))


            C.append('np < min([n for n in E])')  # E is an equivalence class
            B = []  # to replace A

            for f in A:
                if (f.applyMap(np) == np):
                    B.append(f)

            A = B

            # find the largest A-equivalence class

            np = E[0] # first element in the largest equivalence class

        M[n] = C  # M is a list of lists

    return M


#######################################
#### onion decomposition algorithm ####
#######################################

# labels each node with its onion layer and coreness
# G: network to be decomposed
def onionDecompose(G):
    D = sortDegrees(G, G.nodes())
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

            G.remove_node(v)  # delete from D

        layer = layer + 1

        if (D['degree'][0] > core + 1):
            core = D['degree'][0]

    return labels
