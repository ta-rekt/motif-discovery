# this is to try an implementation of a dynamic version of isomorphicExtensions
# Tarek Tohme, July 2019

import numpy as np
import pandas as pd
import networkx as nx
import collections as co
import operator

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


        # equality for functions
    def isEqual(self, g):
        a = set([tuple(k) for k in self.map])
        b = set([tuple(k) for k in g.getMap()])

        return (a == b)


    # finds the inverse of a function
    def inverse(self):
        g = Map([tuple(k)[::-1] for k in self.map])

        return g

    # tests membership of f in a given set
    # speed improvement: try hashing functinons by changing them to strings (on condition that
    # functions in set are sorted by domain)
    def inSet(self, set):
        for f in set:
            if(self.isEqual(f)):
                return True

        return False


    # correctness for bijective functions. checks if function is equal
    # to its inverse
    def isBijection(self):
        g = self.inverse()

        return self.isEqual(g)


    # returns function domain
    def getDomain(self):
        return self.map['domainNode']


    # returns function range
    def getRange(self):
        return self.map['rangeNode']


    # returns domain and range in structured array
    def getMap(self):
        return self.map


    # returns the image of a given argument, -1 if it's not part of the map
    def applyMap(self, node):
        for i in range(len(self.map['domainNode'])):
            if (self.map['domainNode'][i] == node):
                return self.map['rangeNode'][i]
        return -1


class DynamicTree():

    def __init__(self, G, root):
        self.dict = {(root,): {'parent': None,
                            'children': set(),
                            'table': {0: set(G.nodes()).difference(G[root])}},
                    'root': (root,),
                    'leaves': {(root,)}}

        self.graph = G


    def __getitem__(self, key):
        return self.dict[key]


    def __setitem__(self, key, value):
        self.dict[key] = value


    def addNode(self, node, parent):
        if (type(parent) == type(0)):
            parent = (parent,)

        if (node in self.dict.keys()):
            print('node already added')

        else:
            v = node[len(node)-1]
            self.dict[node] = {'parent': parent,
                               'children': [],
                               'table': {0: self.dict[parent]['table'][0].difference(self.graph[v])}}

            self.dict[parent]['children'].add(node)
            self.dict['leaves'].remove(parent)
            self.dict['leaves'].add(node)


    def removeNode(self, node):
        if (node not in self.dict.keys()):
            print('node not in tree')

        else:
            p = self.dict[node]['parent']
            self.dict[p]['children'].remove(node)

            if self.dict[node]['children']:
                for child in self.dict[node]['children']:
                    self.removeNode[child]

            del self.dict[node]

            if node in self.dict['leaves']:
                self.dict['leaves'].remove(node)


    def depth(self, node):
        r = self.dict['root']
        n = node
        path = []

        while(n != r):
            path.append(n)
            n = self.dict[n]['parent']

        return len(path)


    def getTable(self, node):
        return self.dict[node]['table']


    def setTable(self, node, k, value):
        self.dict[node]['table'][k] = value


    def leaves(self):
        return self.dict['leaves']


    def parent(self, node):
        return self.dict[node]['parent']


    def children(self, node):
        return self.dict[node]['children']


    def printTree(self):
        print(self.dict)


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

    if(G.degree(g) >= H.degree(h)):
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



###############################################
#### grochow-kellis motif search algorithm ####
###############################################

def print_situation(D, R, m, l, nmD, fnmD, candidateTree):
    if (len(l) == 1):
        print('D: ', D)

# adapted to iterative dynamic implementation of isomorphicExtensions
def findSubgraphInstances_dynamic(H, G, withSBC=True):

    instances = []

    # adds onion layer, coreness and traversal attributes to G and H
    traversal_order_G = onionDecompose(G)
    traversal_order_H = onionDecompose(H)

    # last node to be peeled off in H
    k = traversal_order_H[len(H)-1]

    # traverse G in increasing onion layer
    for i in range(len(G)):
        g = traversal_order_G[i]

        if(canSupport(k, g, H, G)):
            instances.append(isomorphicExtensions_dynamic(g, H, G))

    return instances


# a dynamic implementation of isomorphic extensions
def isomorphicExtensions_dynamic(g, H, G, M=None):

    instances = []
    traversal = {}

    for i in H.nodes():
        traversal[i] = H.nodes[i]['traversal']

                                 # variables initialization #

    t = len(H)-1                                                     # initial traversal index
    h = traversal[t]                               # last node peeled off is first to traverse
    f = Map([(h, g)])                                                    # initial partial map
    D = {h}                                                           # initial partial domain
    R = {f.applyMap(h)}                                                # initial partial range
    m = h                                                           # initial domain extension
    nmD = D.intersection(set(H[m]))                                      # neighbors of m in D
    fnmD = set([f.applyMap(i) for i in nmD])                          # images of nmD elements
    k = len(nmD)
    c = H.nodes[m]['coreness']                                  # coreness of domain extension

    candidateTree = DynamicTree(G, g)
    leaves = list(candidateTree['leaves'])[:]


                             # computing isomorphism candidates #

    for l in leaves:                              # leaves are in fact paths from root to leaf

        print_situation(D, R, m, l, nmD, fnmD, candidateTree)

        if (type(l) == type(1)):
            v = l
        else:
            v = l[len(l)-1]                                  # actual leaf, candidate for f(m)
        p = candidateTree[l]['parent']                                # parent of current leaf

        if (l == candidateTree['root']):
            candidateTree[l]['table'][1] = set(G[v])
            candidateTree[l]['table'][2] = set()

        else:
            loopEnd = min(c, len(D))+1

            for i in range(1, loopEnd):                                 # [0] is non-neighbors
                N_ip1_R = set(candidateTree[p]['table'][i])
                N_i_R = set(candidateTree[p]['table'][i-1])
                N_v = set(G[v])

                newEntry = N_ip1_R.union(N_i_R.intersection(N_v)).difference(N_ip1_R.intersection(N_v))
                candidateTree[l]['table'][i] = newEntry

            candidateTree[l]['table'][loopEnd] = set()

        N_k = candidateTree[l]['table'][k]                              # get N_k neighborhood


                                   # testing isomorphism #

        if not N_k:                                                # N_k neighborhood is empty
            if not (l == candidateTree['root']):
                while (len(candidateTree[p]['children']) == 1):
                    f = p
                    p = candidateTree[p]['parent']

                candidateTree.removeNode(f)

        else:
            print(N_k)
            for i in N_k:
                if (set(G[i]).intersection(D) == fnmD):                            # O(|G[i]|)
                    temp = list(l)
                    temp.append(i)
                    newNode = tuple(temp)
                    candidateTree.addNode(newNode, l)

        if (l == leaves[len(leaves)-1]):                                        # last element
            if (len(D) == len(H)):                                                  # all done
                for leaf in leaves:
                    instances.append[leaf]

                return instances

            D = D.union({m})
            t = t - 1
            m = traversal[t]
            nmD = D.intersection(set(H[m]))                              # neighbors of m in D
            fnmD = set([f.applyMap(i) for i in nmD])                  # images of nmD elements
            k = len(nmD)
            c = H.nodes[m]['coreness']
            leaves = list(candidateTree['leaves'])[:]


#######################################
#### onion decomposition algorithm ####
#######################################

# labels each node with its onion layer and coreness, adds them as attributes to G
# G: network to be decomposed
def onionDecompose(G):
    K = nx.Graph()
    K.add_edges_from(G.edges())

    coreness = {}
    onion_layer = {}
    traversal_order = {}

    core = 1
    layer = 1
    count = 0

    while (len(K.nodes()) > 0):
        thisLayer = [v for v in K.nodes() if len(K[v]) <= core]

        for v in thisLayer:
            coreness[v] = core
            onion_layer[v] = layer
            traversal_order[count] = v

            K.remove_node(v)  # delete from D
            count = count + 1

        layer = layer + 1

        D = [len(K[k]) for k in K.nodes()]

        if(D):
            minn = min(D)

            if (minn > core):
                core = minn
    nx.set_node_attributes(G, coreness, 'coreness')
    nx.set_node_attributes(G, onion_layer, 'onion_layer')
    nx.set_node_attributes(G, traversal_order, 'traversal')

    return traversal_order





























# finds symmetry-breaking conditions for H given HE and Aut(H). don't need a labeling function on G
# because nodes are already numbered. just need to find the conditions.
# you wouldn't need an object of type condition: a set of nodes is enough.
# given a node n of H and the equivalence class of n as determined by symmetryConditions (not
# the actual equivalence class) one only has to check if the newly encountered node in G has a
# label less than all the labels of the images of the nodes in that equivalence class, which would
# be stored alongside n in M. if a node has no equivalence class in M, then there are no symmetry-
# breaking constraints on that node.
# HE: set of equivalence representatives of H. these are nodes of H, one from each equivalence class.
# aut: set of automorphismisms of H. each automorphism is an object of type Map that has identical
# domain and range, and where all nodes are distinct and appear exactly once
def symmetryConditions(aut):
    M = {} # dict containing nodes in H that have conditions and the set of nodes that constrain them
    eqClassesAndReps = findEquivalenceClasses(aut) # list of tuples: (eq. class, representative node)
    HE = [t[1] for t in eqClassesAndReps]                         # list of representative nodes only

    for i in range(len(HE)):

        n = HE[i]
        np = n
        A = aut
        S = eqClassesAndReps[i][0].difference({n})                   # equivalence class of n minus n
        M[n] = S   # S is a set of nodes such that l(n) < min(l(k) | k in S) i.e. l(S) must be > l(n)

        while len(A) > 1:

            A = [f for f in A if f.applyMap(np) == np]       # "pinch" the set of automorphisms at np

            tempEqClassesAndReps = findEquivalenceClasses(A)

            a = [len(t[0]) for t in tempEqClassesAndReps]
            maxx = max(a)                                      # max size of equivalence classes in A

            ind = [key for key, val in enumerate(a) if val == maxx]
            ind = ind[0]           # argmax. np.argmax() was being a pain in the neck for some reason

            np = tempEqClassesAndReps[ind][1]           # representative of largest equivalence class
            Sp = tempEqClassesAndReps[ind][0].difference({np})
            M[np] = Sp

    return M


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


# returns only functions that are unique in the set maps
# maps: set of maps, possibly with duplicates
def bijectionsOnly(maps):
    # dummy = []
    to_return = []

    for f in maps:
        if (f.isBijection()):
            to_return.append(f)

    return to_return


# there's 4 possible cases:
# 1. m is a key: we should check if all the images of the values of m
#    in M are greater than n. if so, accept n. if not, reject the
#    whole instance (figure out how to do that). if some values have
#    not been mapped yet, accept n.
# 2. m is a value: we should find its corresponding key and check if
#    its image is less than n.
# 3. m is a key, but has no values: accept n.
# 4. m is neither a key nor a value: accept n.
# m: domain extension node
# n: range extension node
# M: symmetry-breaking conditions
# f: partial map
def checkSBC(m, n, M, f):
    if (M == None):
        return True

    if (m in M.keys()):
        value = M[m]

        if (value):  # case 1
            images = [f.applyMap(k) for k in value] # {L(k)|k in values}

            if (n < min(images)): # all values mapped
                # printSBC(m, n, M, f)
                # print('min(images): ',end='')
                # print(min(images))

                if (m in k for k in M.values()):   # case 2
                    keys = [list(M.keys())[list(M.values()).index(k)] for k in M.values() if m in k]

                    for k in keys:
                        if (f.applyMap(k) < 0):
                            # print('some conditions unknown. node accepted')
                            # printSBC(m, n, M, f)
                            return True

                        else:
                            if (f.applyMap(k) > n):
                                # print('violates to sbc. node rejected')
                                # printSBC(m, n, M, f)
                                # print('f(key of m)')
                                # print(f.applyMap(key))
                                return False

                    # print('conforms to sbc. node accepted')
                    # printSBC(m, n, M, f)
                    # print('f(key of m)')
                    # print(f.applyMap(keys))
                    # print('min(images): ',end='')
                    # print(min(images))

                    return True


            elif(min(images) == -1): # some values not mapped
                # print('some conditions unknown. node accepted')
                # printSBC(m, n, M, f)
                return True

            else:  # n violates sbc
                # print('violates sbc. node rejected')
                return False


    if (m in k for k in M.values()):   # case 2
        keys = [list(M.keys())[list(M.values()).index(k)] for k in M.values() if m in k]

        for k in keys:
            if (f.applyMap(k) < 0):
                # print('some conditions unknown. node accepted')
                # printSBC(m, n, M, f)
                return True

            else:
                if (f.applyMap(k) > n):
                    # print('violates to sbc. node rejected')
                    # printSBC(m, n, M, f)
                    # print('f(key of m)')
                    # print(f.applyMap(key))
                    return False

        # print('conforms to sbc. node accepted')
        # printSBC(m, n, M, f)
        # print('f(key of m)')
        # print(f.applyMap(keys))

        return True

    else:  # case 4
        # print('no conditions on n. node accepted')
        # printSBC(m, n, M, f)
        return True

    # print('something wrong with SBC')
    # printSBC(m, n, M, f)
    return False


# to ease debugging
def printSBC(m, n, M, f):
    print('m: ',end='')
    print(m)
    print('n: ',end='')
    print(n)
    print('M: ',end='')
    print(M)
    print('f: ',end='')
    print(tuple(f.getMap()))
