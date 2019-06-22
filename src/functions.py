# this is for testing implementations of graph algorithms to speed up motif seach in a network.
# after refactoring, testing, verification, and new additions make this available as open-source.
# also, the comments have to be organized using a documentation tool. i like to throw comments
# everywhere becaue it gives consolation that i'm writing useful stuff.
# Tarek Tohme, June 2019

import numpy as np
import pandas as pd
import networkx as nx
import collections as co

count = 0

######################################
#### useful functions and classes ####
######################################

# function f: V -> V that maps nodes in H to nodes in G
# f: ordered list of tuples representing the partial map between D and R, respectively the domain and
# range of f. the function makes sure it is a bijection each time an extension is attempted. if a new
# node is the same as some existing node, the function replaces the old node with the new node.
class Map():

    def __init__(self, init):   # init is a list containing the initial nodes of the partial map
        self.map = np.array(init, dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])

    def extend(self, extension):   # extends the partial map by a list of tuples called extension
        # print('##### EXTENDING PARTIAL MAP #####')
        # print('extension: ', end='')
        # print(extension)

        temp = np.array(extension, dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])
        # print(self.map.shape)
        # print(temp.shape)
        map = self.map

        if (not extension):
            # print('extension is null')
            return False

        # checking for duplicates within extension
        for i in temp:
            l = [m for m,n in enumerate(temp) if tuple(n)==tuple(i)]

            if (len(l) > 1):
                # print('cannot extend function. duplicate element in extension')
                self.map = map
                return False

            for j in temp:

                if(np.array_equal(i, j) == False):

                    if (i['domainNode'] == j['domainNode']):
                        # print('cannot extend function. duplicate domain node in extension')
                        self.map = map
                        return False

                    if (i['rangeNode'] == j['rangeNode']):
                        # print('cannot extend function. duplicate range node in extension')
                        self.map = map
                        return False

            # checking for duplicates between extension and partial map
            for j in self.map:

                if (i['domainNode'] == j['domainNode']):
                    # print('duplicate domain node. using new node.')
                    duplicate = i['domainNode']
                    b = [x for k,x in enumerate(self.map) if self.map['domainNode'][k] != duplicate]
                    self.map = np.array(b, dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])

                elif (i['rangeNode'] == j['rangeNode']):
                    # print('duplicate range node. using new node.')
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


# returns true if g can support h (according to node degree and neighbor degree sequence)
# and false otherwise.
# g, h: indices of corresponding nodes
# H: query graph
# G: network to be queried
def canSupport(h, g, H, G):
    # print('##### CAN SUPPORT CALL #####')
    # print('query node: ', end='')
    # print(h)
    # print('network node: ', end='')
    # print(g)

    if(G.degree[g] >= H.degree(h)):
        neighborsH = sortDegrees(H, H[h])
        neighborsG = sortDegrees(G, G[g])
        neighborsH = np.array(list(reversed(neighborsH)), dtype=[('node', 'i4'), ('degree', 'i4')])
        neighborsG = np.array(list(reversed(neighborsG)), dtype=[('node', 'i4'), ('degree', 'i4')])

        minn = min(len(neighborsH), len(neighborsG))

        for i in range(minn):
            if(neighborsG['degree'][i] < neighborsH['degree'][i]):
                # print("rejected because of sequence")
                return False
        # print("can support")
        return True

    # print("rejected because of degree")
    return False


# returns an ordered list of tuples (node, degree) sorted by degree
# H: network
# L: a list of nodes that should be sorted
def sortDegrees(H, L):
    # print('##### SORT DEGREES CALL #####')
    sortedDegreeH = H.degree()
    sortedDegreeH = [(i, sortedDegreeH(i)) for i in L]
    sortedDegreeH = np.array(sortedDegreeH, dtype=[('node', 'i4'), ('degree', 'i4')])
    sortedDegreeH = np.sort(sortedDegreeH, order='degree')
    # print('sorted degrees: ', end='')
    # print(sortedDegreeH)
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
    # print('##### MOST CONSTRAINED NODE #####')
    # print('domain of partial map: ', end='')
    # print(D)

    maxx = 0
    candidates = np.array([], dtype=[('node', 'i4'), ('attribute', 'i4')])
    neighborsD = {None}

    # list of neighbors of D
    for d in D:
        neighborsD = neighborsD.union(set(H[d]))

    neighborsD = neighborsD.difference({None})

    # neighbors of D exclusively in H with the number of neighbors in D
    neighborsD = neighborsD.difference(set(D))
    # print('neighbors of D: ', end='')
    # print(neighborsD)

    for n in neighborsD:
        count = 0
        for i in H[n]:
            if (i in D):
                count = count + 1

        c = np.array([(n, count)], dtype=[('node', 'i4'), ('attribute', 'i4')])
        candidates = np.concatenate((candidates, c))

    # select the candidates with the most neighbors in D
    # print('candidates: ', end='')
    # print(candidates)

    maxx = max(candidates['attribute'])
    indices = [key for key, val in candidates if val == maxx]

    # degrees of nodes with most neighbors in D
    candidates = np.array([(i, len(H[i])) for i in indices], dtype=[('node', 'i4'), ('attribute', 'i4')])

    if (len(candidates) == 1):
        return candidates['node'][0]

    else:
        # select those with highest degree
        # print('nodes with most neighbors in D, degrees: ',end='')
        # print(candidates)

        maxx = max(candidates['attribute'])
        ind = [key for key, val in candidates if val == maxx]
        # print('nodes of largest degree: ', end='')
        # print(ind)

        if (len(ind) == 1):
            return ind[0]

        else:
            # select those with highest degree sequence
            m = largestDegreeSequence(H, ind)
            # print('most constrained node: ', end='')
            # print(m)
            return m[0]


# finds the set of equivalence classes. returns a list; the position
# of an element in the list corresponds to a node and that element is
# the equivalence class it belongs to.
# aut: set of automorphisms. must be a list of objects of type Map.
def findEquivalenceClasses(aut):
    eq = []

    for n in aut[0].getDomain():
        E = {None}

        for f in aut:
            E = E.union({f.applyMap(n)})

        E = E.difference({None})

        tuple = (E, list(E)[0])
        if (tuple not in eq):
            eq.append(tuple)

    return eq


###############################################
#### grochow-kellis motif search algorithm ####
###############################################

# finds all instances of query graph H in network G
# H: query graph
# G: network to be queried
def findSubgraphInstances(H, G):
    # print('##### FIND SUBGRAPH INSTANCES #####')

    # instances of H found in G
    instances = []

    # sort nodes of G by degree
    sortedDegreeG = sortDegrees(G, G.nodes())

    # sort by degree sequence. larger sequence is pointwise larger for all nodes (later)

    for g in sortedDegreeG['node']:
        for h in H.nodes():
            if(canSupport(h, g, H, G)):
                f = Map([(h, g)])  # initialize partial map with f(h) = g. argument must be list
                iso = isomorphicExtensions(f, H, G)

                if(iso):  # sometimes iso is empty
                    instances.append(iso)

        # G.remove_node(g)
    instances = [next(iter(k)) for k in instances]  # removes the elements from their sets

    return instances

# finds all isomorphic extensions of a partial map [satisfying the symmetry-breaking
# condition C at h]. returns them in a list of tuples [(a1, b1), ... ,(ak, bk)]
# f: partial map to be extended
def isomorphicExtensions(f, H, G): #, C, h)
    print('##### ISOMORPHIC EXTENSIONS CALL #####')


    isomorphisms = set()
    neighborsR = {None}
    neighborsD = {None}
    D = f.getDomain()
    R = f.getRange()

    if(set(D) == set(H.nodes())):
        print('INSTANCE FOUND')
        print(list(f.getMap()))
        f.extend(list(np.sort(f.getMap(), order='domainNode')))
        return {f}

    m = mostConstrainedNode(D, H)
    print('domain extension (m): ', end='')
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
    # print('neighbors of partial domain: ',end='')
    # print(neighborsD)

    # check for induced isomorphism
    for n in neighborsR:

        print('partial map: ')
        print(np.vstack(f.getMap()))

        print('try range extension (n): ',end='')
        print(n)

        neighbMinD = set(H[m]).intersection(set(f.getDomain()))
        f_neighbMinD = set([f.applyMap(k) for k in neighbMinD])
        neighbNinR = set(G[n]).intersection(set(f.getRange()))

        # print('neighbors of m in D: ', end='')
        # print(neighbMinD)
        # print('neighbors of n in R: ', end='')
        # print(neighbNinR)
        # print('f(neighbors of m in D): ', end='')
        # print(f_neighbMinD)

        if(set(f_neighbMinD) != set(neighbNinR)):
            print('range extension not valid')
            pass

        else:
            print('range extension valid')
            fp = Map(list(f.getMap()))
            newNode = np.array((m, n), dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])

            fp.extend([newNode])
            print('called')
            iso = isomorphicExtensions(fp, H, G)
            print('returned')
            isomorphisms = isomorphisms.union(iso)

            # if(isomorphisms):
            #     print(next(iter(isomorphisms)).getMap())

    # print('EXIT ISOMORPHIC EXTENSIONS')
    # print('isomorphicExtensions output: ',end='')
    # print(isomorphisms)

    return isomorphisms


# finds symmetry-breaking conditions for H given HE and Aut(H). don't need a labeling function on G
# because nodes are already numbered. just need to find the conditions.
# HE: set of equivalence representatives of H. these are nodes of H, one from each equivalence class.
# aut: set of automorphismisms of H. each automorphism is an object of type Map that has identical
# domain and range, and where all nodes are distinct and appear exactly once

# additional comments:
# you wouldn't need an object of type condition: a set of nodes is enough.
# given a node n of H and the equivalence class of n as determined by symmetryConditions (not
# the actual equivalence class) one only has to check if the newly encountered node in G has a
# label less than all the labels of the images of the nodes in that equivalence class, which would
# be stored alongside n in M. if a node has no equivalence class in M, then there are no symmetry-
# breaking constraints on that node.
def symmetryConditions(aut):
    M = {}  # dict containing nodes in H that have conditions and the set of nodes that constrain them
    eqClassesAndReps = findEquivalenceClasses(aut)  # list of tuples: (eq. class, representative node)
    HE = [t[1] for t in eqClassesAndReps]  # list of representative nodes only

    for i in range(len(HE)):

        n = HE[i]
        np = n
        A = aut
        S = eqClassesAndReps[i][0].difference({n})  # equivalence class of n minus n
        M[n] = S  # S is a set of nodes such that l(n) < min(l(k) | k in S) i.e. l(S) must be > l(n)

        while len(A) > 1:

            A = [f for f in A if f.applyMap(np) == np]  # "pinch" the set of automorphisms at np

            tempEqClassesAndReps = findEquivalenceClasses(A)

            a = [len(t[0]) for t in tempEqClassesAndReps]
            maxx = max(a)  # max size of equivalence classes in A

            ind = [key for key, val in enumerate(a) if val == maxx]
            ind = ind[0]   # argmax. np.argmax() was being a pain in the neck for some reason

            np = tempEqClassesAndReps[ind][1]  # representative of largest equivalence class
            Sp = tempEqClassesAndReps[ind][0].difference({np})
            M[np] = Sp

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
