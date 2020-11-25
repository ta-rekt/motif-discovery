# this is for testing implementations of graph algorithms to speed up motif seach in a network.
# after refactoring, testing, verification, and new additions, the intent is to contribute this to networkx.
# also, the comments have to be organized using a documentation tool.
# Tarek Tohme, June 2019

import numpy as np
import pandas as pd
import networkx as nx
import collections as co
import random as rnd
import permanent_fast as per_f
import operator
from operator import mul
import itertools as it
from itertools import groupby
import functions_regular_wc as fnr_wc
import math
import scipy.special as sp

count = 0

##########################################################################
###################### useful functions and classes ######################
##########################################################################

# function f: V -> V that maps nodes in H to nodes in G
# f: ordered list of tuples representing the partial map between D and R, respectively the domain and
# range of f. the function makes sure it is a bijection each time an extension is attempted. if a new
# node is the same as some existing node, the function replaces the old node with the new node.
class Map():

    def __init__(self, init, multiplier=1):   # init is a list containing the initial nodes of the partial map
        self.map = np.array(init, dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])
        self.multiplier = multiplier

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

    # sets the multiplier
    def setMult(self, mult):
        self.multiplier = mult

    # gets the multiplier
    def getMult(self):
        return self.multiplier


# returns true if g can support h (according to node degree and neighbor degree sequence)
# and false otherwise.
# g, h: indices of corresponding nodes
# H: query graph
# G: network to be queried
def canSupport(h, g, H, G):
    # print('##### CAN SUPPORT CALL #####')
    # print('query node: ', end='')
    # print(h, H[h])
    # print('network node: ', end='')
    # print(g, G[g])

    if (len(G[g]) >= len(H[h])):
        neighborsH = sortDegrees(H, H[h])
        neighborsG = sortDegrees(G, G[g])
        neighborsH = np.array(list(reversed(neighborsH)), dtype=[('node', 'i4'), ('degree', 'i4')])
        neighborsG = np.array(list(reversed(neighborsG)), dtype=[('node', 'i4'), ('degree', 'i4')])

        minn = min(len(neighborsH), len(neighborsG))

        for i in range(minn):
            if (neighborsG['degree'][i] < neighborsH['degree'][i]):
                # print("rejected because of sequence")
                return False

        if (H.node[h]['coreness'] > G.node[g]['coreness']):
                # print("rejected because of coreness")
                return False

        # print("can support")
        return True

    # print("rejected because of degree")
    # print(h, g)
    # print('N(h)', list(H[h]), 'N(g)', list(G[g]))
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
    candidates = np.array([(i, len(H[i])) for i in indices], dtype=[('node','i4'),('attribute','i4')])

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


# finds the set of equivalence classes. returns a dict; keys are representative nodes
# and values are equivalence classes
# aut: set of automorphisms. must be a list of objects of type Map.
def findEquivalenceClasses(aut, M=None):
    eq = {}

    for m in aut[0].getDomain():   # loop through all nodes in the domain (in order)
        E = set()

        for f in aut:
            n = f.applyMap(m)

            if (checkSBC(m, n, M, f)):
                E = E.union({n})   # gather all nodes that n can be mapped to in E

        if (E not in eq.values()):
            eq[m] = E

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


# finds the neighbors of m in D, their images, and the neighbors of n in R
# comments: using the dynamic algorithm, nmd is in O(l), the coreness of the shell of m in H,
# and nnr is computed dynamically in O(l)
def findCandidates(m, n, H, G, f):
    nmd = set(H[m]).intersection(set(f.getDomain())) # neighbors of m in D
    fnmd = set([f.applyMap(k) for k in nmd])
    nnr = set(G[n]).intersection(set(f.getRange()))  # neighbors of n in R

    return [nmd, fnmd, nnr]


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
                                # print(f.applyMap(k))

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

# relabels a graph's edges in increasing degree order
# H: a networkx graph
def relabelByDegree(H):
    sortedDegreeH = sortDegrees(H, H.nodes())
    mappingH = dict(zip(sortedDegreeH['node'], range(len(H))))
    H = nx.relabel_nodes(H, mappingH)

    return H


###################################################################################
###################### grochow-kellis motif search algorithm ######################
###################################################################################

def findSubgraphInstances(strH, strG):
    instances = []
    shell_1_counts = 0

    H = nx.Graph()
    G = nx.Graph()
    H = parseStringToGraph(strH)
    G = parseStringToGraph(strG)

    H.remove_edges_from(nx.selfloop_edges(H))
    G.remove_edges_from(nx.selfloop_edges(G))

    G_order = onionDecompose(G)
    onionDecompose(H)

    marks_G = {}
    marks_H = {}

    K = nx.Graph()
    K.add_edges_from(H.edges())

    attributesH = H.nodes[0].keys()
    for name in attributesH:
        attr = nx.get_node_attributes(H, name)
        nx.set_node_attributes(K, attr, name)

    aut = fnr_wc.findSubgraphInstances(H, K, False)
    eqClasses = fnr_wc.findEquivalenceClasses(aut)
    HE = list(eqClasses.keys())
    M = fnr_wc.symmetryConditions(aut)
#     print('M', M)

    is_tree = nx.is_tree(H)

    combos = {}

    for h in HE:
        if (is_tree):
            combos[h] = enumSubtrees(H, h, True)

        elif (H.nodes[h]['coreness'] == 1):
            p = -1

            for i in H[h]:
                if (H.nodes[i]['onion_layer'] > H.nodes[h]['onion_layer']):
                    p = i

            T = subtree(H, h, p)
            combos[h] = enumSubtrees(T, h, True)

        elif (H.nodes[h]['is_root']):
            T = detachTree(H, h)
            T = T[0]
            combos[h] = enumSubtrees(T, h, True)


    for c in range(len(G)):
        g = G_order[c]

        if (g in G):
            # print('g:', g, 'coreness:', G.nodes[g]['coreness'])

            if (G.nodes[g]['coreness'] == 1):
                # print('1-SHELL')

                if (is_tree):
                    p = -1

                    for i in G[g]:
                        if (G.nodes[i]['onion_layer'] > G.nodes[g]['onion_layer']):
                            p = i

                    if (p >= 0):

                        S = subtree(G, g, p)

                        for h in HE:
                            count = countRootedSubtrees(H, h, S, g)

                            if (count > 0):
                                shell_1_counts += count
                                # print('root:', g, 'count:', shell_1_counts)


            elif (G.nodes[g]['is_root']):
                # print('IN-BETWEEN')

                S = detachTree(G, g)
                S = S[0]
                # print('S:', S.nodes())

                for h in HE:
                    if (H.nodes[h]['coreness'] == 1 or H.nodes[h]['is_root'] or is_tree):
                        # print('h:', h)
                        combo = combos[h]
                        # print('combo:', combo)

                        for i in combo.values():
                            for select in i:
                                # print('select:', select)

                                if (select):
                                    T = detachTree(H, h, select)
                                    A = T[0]
                                    B = T[1]
                                    # print('A:', A.nodes())

                                    f = fnr_wc.Map([(h, g)])

                                    if (len(S) >= len(A)):

                                        count = countRootedSubtrees(A, h, S, g)
                                        # print('count:', count)
                                        f.setMult(count)

                                        if (count > 0):

                                            for s in S:
                                                if (s != g):
                                                    marks_G[s] = True
    #                                                 print(s, 'marked')

                                            nx.set_node_attributes(G, marks_G, 'marked')

                                            # print('f before isoext:', f.getMap())
                                            N = copy.deepcopy(M)

                                            for key in M.keys():
                                                for val in M[key]:
                                                    if (val not in B):
                                                        N[key].remove(val)

                                                if (key not in B):
                                                    del N[key]

                                            iso = fnr_wc.isomorphicExtensions(f, B, G, 0, N)

                                            for s in S:
                                                if (s != g):
                                                    marks_G[s] = False
    #                                                 print(s, 'unmarked')

                                            nx.set_node_attributes(G, marks_G, 'marked')

                                            if(type(iso) == type(f)):
                                                instances.append(iso)

                                            else:
                                                for i in iso:
                                                    instances.append(i)

                                            # print('instances:', [[tuple(x.getMap()), x.getMult()] for x in instances])

                for s in S.nodes():
                    if (G.nodes[s]['coreness'] == 1):
                        G.remove_node(s)

    # for g in list(G.nodes()):
    #     if (G.nodes[g]['coreness'] == 1):
    #         G.remove_node(g)

    # print('1-shell:', [[tuple(f.getMap()), f.getMult()] for f in instances])

    total = 0
    for f in instances:
        total = total + f.getMult()

    # print('2-CORE')
    two_core = fnr_wc.findSubgraphInstances(H, G, True)

    # print('1-shell:', shell_1_counts)
    # print('in-between:', total)
    # print('2-core:', two_core)

    return int(total + two_core + shell_1_counts)


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
# M[n] = S: S is a set of nodes such that l(n) < min(l(k) | k in S) i.e. l(S) must be > l(n)
def symmetryConditions(aut):
    M = {}  # dict containing nodes in H that have conditions and their constraining sets of nodes
    eqClasses = findEquivalenceClasses(aut)  # list of tuples: eq. class, representative node
    HE = [i for i in eqClasses.keys()]  # list of representative nodes

    # print(eqClasses)

    for n in HE:
        if (n not in M.keys()):
            np = n
            A = aut
            temp_eqClasses = eqClasses

            # print(n)

            A = [f for f in A if checkSBCFunction(f, M)]

            while (len(A) > 1):

                # print(len(A))
                # print('A before:')
                # print([tuple(i.getMap()) for i in A])

                Sp = temp_eqClasses[np]
                Sp.remove(np)
                M[np] = Sp

                # print('M[',np,'] = ',Sp)

                A = [f for f in A if f.applyMap(np) == np and checkSBCFunction(f, M)]
                temp_eqClasses = findEquivalenceClasses(A, M)

                # print(temp_eqClasses)

                maxx = 0

                for k in temp_eqClasses.keys():
                    size = len(temp_eqClasses[k])

                    if (size > maxx):
                        maxx = size
                        np = k

                # print('A after:')
                # print([tuple(i.getMap()) for i in A])

    return M


# checks if a function satisfies a set of SBCs
def checkSBCFunction(f, M):
    for m in f.getDomain():
        if (not checkSBC(m, f.applyMap(m), M, f)):
            return False

    return True



###########################################################################
###################### onion decomposition algorithm ######################
###########################################################################

# labels each node with its onion layer and coreness, adds them as attributes to G
# G: network to be decomposed
def onionDecompose(G):
    K = nx.Graph()
    K.add_edges_from(G.edges())

    is_root = {}
    marked = {}
    coreness = {}
    onion_layer = {}
    traversal_order = {}

    core = 1
    layer = 1
    count = 0

    shell_2 = set()

    while (len(K.nodes()) > 0):
        thisLayer = [v for v in K.nodes() if len(K[v]) <= core]

        for v in thisLayer:
            coreness[v] = core
            onion_layer[v] = layer
            traversal_order[count] = v
            is_root[v] = []

            if (core >= 2):
                shell_2.add(v)

            K.remove_node(v)  # delete from D
            count = count + 1

        layer = layer + 1

        D = [len(K[k]) for k in K.nodes()]

        if(D):
            minn = min(D)

            if (minn > core):
                core = minn

    for v in G:
        marked[v] = False

        if (v in shell_2):
            for u in G[v]:
                if coreness[u] == 1:
                    is_root[v].append(u)

    nx.set_node_attributes(G, marked, 'marked')
    nx.set_node_attributes(G, is_root, 'is_root')
    nx.set_node_attributes(G, coreness, 'coreness')
    nx.set_node_attributes(G, onion_layer, 'onion_layer')
    nx.set_node_attributes(G, traversal_order, 'traversal')

    return traversal_order




##########################################################################
###################### counting isomorphic subtrees ######################
##########################################################################

# generates a random tree of size n
def randomTree(n, seed):
    rnd.seed(seed)

    S = nx.Graph()
    S.add_node(0)

    for i in range(1, n):
        r = rnd.randint(0,i-1)
        S.add_edge(r, i)

    sortedDegreeS = sortDegrees(S, S.nodes())
    mapping = dict(zip(sortedDegreeS['node'], range(len(S))))
    S = nx.relabel_nodes(S, mapping)

    return S


# returns the preorder traversal sequence of a tree. for initial call, parent = root
def compressTree(T, node, parent, dict = {}):

    if (not nx.is_tree(T)):
        return False

    string = ''
    root = -1

    if(node == parent):
        children = set(T[parent])
        root = node

    else:
        children = set(T[node]).difference({parent})

    if(children):
        for n in children:
            s = '1'
            substring = compressTree(T, n, node, dict)
            dict[n] = substring
            s = s + substring
            s = s + '0'
            string = string + s

    else:
        return ''

    if (root > -1):
        dict[root] = string
        return dict

    return string


# returns a list of depths with the sorted degree sequence of the
# vertices at that depth
def depthDegreeSequence(T, root):
    paths = nx.single_source_shortest_path(T, root)
    seq = {}
    longestPath = max([len(v) for v in paths.values()])

    for i in range(1, longestPath+1):
        seq[i] = set()

    for key in paths.keys():
        path = paths[key]
        l = len(path)
        seq[l] = seq[l].union({key})

    for key in seq.keys():
        degrees = sorted([len(T[v])-1 for v in seq[key] if key != 1])
        seq[key] = [i for i in reversed(degrees)]

    seq[1] = [len(T[root])]

    return seq


# counts the matchings between every subtree of T with S and stores
# them in the dictionary
def countSubtreeMatchings(T, rootT, parentT, S, rootS, dic):
    dic[rootT] = countRootedSubtrees(T, rootT, S, rootS)

    if (rootT == parentT):
        children = set(T[rootT])
    else:
        children = set(T[rootT]).difference({parentT})

    if (children):
        for i in children:
            Ti = detachTree(T, i, rootT)
            countSubtreeMatchings(Ti, i, rootT, S, rootS, dic)
    return


# returns the subtree of T rooted at r with parent p
# can work on non-trees, as long as F is a tree
def subtree(T, r, p):
    if (r == p):
        return T

    K = nx.Graph()
    K.add_edges_from(T.edges())
    K.remove_edge(r, p)

    for k in nx.connected_components(K):
        if (r in k):
            F = nx.Graph()

            if (len(k) == 1):
                F.add_node(r)

            else:
                for i in k:
                    for j in k:
                        if ((i, j) in K.edges()):
                            F.add_edge(i, j)

    return F


# returns a tree in the 1-shell rooted at root.
# can select which children to include using select.
# use select to properly separate subtrees and 2-core
def detachTree(G, root, select=None):
    K = nx.Graph()
    K.add_edges_from(G.edges())

    children = set()

    for i in G[root]:
        if (G.nodes[i]['coreness'] == 1):
            if (not select or i in select):
                children.add(i)

    T = nx.Graph()
    T.add_node(root)

    for i in children:
        Ti = subtree(K, i, root)
        T.add_edges_from(Ti.edges())
        T.add_edge(root, i)

    for t in T:
        if (t != root):
            K.remove_node(t)

    return [T, K]


# computes the permanent of a non-square matrix
def rectPermanent(M, l, r):
    # obtain the biadjacency matrix B from M, assuming M is
    # properly indexed with left vertices numbered less than
    # right vertices
    ind = [[i+j*(l+r) for i in range(l, l+r)] for j in range(l)]
    B = np.array(np.take(M, ind))

    # enumerate all square submatrices of B
    indices = [i for i in range(r)]
    p = 0

    for e in it.combinations(indices, l):
        ind = [[i+j*(r) for i in e] for j in range(l)]
        A = np.array(np.take(B, ind))
        p += per_f.permanent_fast(A)

    return p


# draws a bipartite graph with labeled edges
def drawBipartiteGraph(F, indices):

    M = nx.to_numpy_matrix(F)
    D = nx.to_dict_of_dicts(F)

    pos = nx.bipartite_layout(F, indices)
    nx.draw_networkx_nodes(F, pos)
    nx.draw_networkx_edges(F, pos, width=1.0, alpha=0.5)
    _dict_ = {}
    _dict2_ = {}

    for n in F.edges():
        if (D[n[0]][n[1]]['weight']):
            _dict_[n] = D[n[0]][n[1]]['weight']

    for n in F.nodes():
        _dict2_[n] = n

    nx.draw_networkx_edge_labels(F, pos, edge_labels=_dict_, font_size=8, label_pos=0.25)
    nx.draw_networkx_labels(F, pos, labels=_dict2_, font_size=8)


# recursively counts the matchings of a bipartite graph G, obtained
# from subtree isomorphisms between subtrees of T and subtrees of S
def countMatchings(T, rootT, S, rootS, preorderTree={}, init=True):
    # print('countmatchings')

    parentT = rootT
    parentS = rootS
    childrenT = set(T[rootT]).difference({rootT})
    childrenS = set(S[rootS]).difference({rootS})
    G = nx.Graph()
    m = 1

    # preorderTree is a dict: keys are roots of each subtree, values are their preorder strings
    if (init):
        preorderTree = compressTree(T, rootT, rootT)

    treeChildrenT = [child for child in childrenT if preorderTree[child] != '']
    starChildrenT = [child for child in childrenT if preorderTree[child] == '']

    treeChildrenS = [child for child in childrenS if len(S[child]) > 1]

    preorderTreeChildren = sorted([preorderTree[i] for i in treeChildrenT])
    classSizes = [len(list(group)) for key, group in groupby(preorderTreeChildren)]

    # print('T', childrenT, treeChildrenT, starChildrenT)
    # print('S', childrenS, treeChildrenS)

    for val in classSizes:
        if (val > 1):
            m *= math.factorial(val)

    # base case
    if (not childrenT):
        return 1

    if (len(treeChildrenT) == 0):
        k_binom = len(childrenT)
        n_binom = len(childrenS)

        return sp.binom(n_binom, k_binom)

    if (len(childrenT) > len(childrenS) or len(treeChildrenT) > len(treeChildrenS)):
        return 0

    indices = {}
    left = set()

    l = len(treeChildrenT)
    r = len(treeChildrenS)

    for key_i, i in enumerate(treeChildrenT):
        for key_j, j in enumerate(treeChildrenS):
            Ti = subtree(T, i, rootT)
            Sj = subtree(S, j, rootS)

            k = countMatchings(Ti, i, Sj, j, preorderTree, False)

            if (k > 0):
                G.add_edge(key_i, key_j+len(treeChildrenT), weight = k)

            indices[key_j+len(treeChildrenT)] = str(j) + 'r'

        indices[key_i] = str(i) + 'l'
        left.add(indices[key_i])

    K = nx.relabel_nodes(G, indices)

    mat = nx.to_numpy_matrix(G, [i for i in range(l+r)])
    matchings = rectPermanent(mat, l, r)

    k_binom = len(starChildrenT)
    n_binom = len(childrenS) - len(treeChildrenT)

    binomCount = sp.binom(n_binom, k_binom)

    return (matchings / m) * binomCount


# recursively find the number of induced subtrees of H's depth-1
# branches and apply bipartite matching to count the number of
# possible matches.
def countRootedSubtrees(T, rootT, S, rootS):
    count = 0

    # early abort by root degree
    if (len(T[rootT]) > len(S[rootS])):

        # print('early abort: root degree')

        return count

    # early abort by depth
    if (nx.eccentricity(T, rootT) > nx.eccentricity(S, rootS)):

        # print('early abort: max depth')

        return count

    seqT = depthDegreeSequence(T, rootT)
    seqS = depthDegreeSequence(S, rootS)

    for i in seqT.keys():
        # early abort by number of nodes at each depth
        if (len(seqT[i]) > len(seqS[i])):

            # print('early abort: number of depth-k nodes')

            return count

        # early abort by node degree sequence at each depth
        for k, v in enumerate(seqT[i]):
            if (v > seqS[i][k]):

                # print('early abort: depth-k degree sequence')

                return count

    # preorder string of T
    preorder = compressTree(T, rootT, rootT)

    # recursion
    count = countMatchings(T, rootT, S, rootS, preorder)

    return count


# enumerates all combinations of h's children in T up to symmetry
def enumSubtrees(T, h, sb=False):
    combos = {}

    children = list(T[h])

#     print('T', T.nodes())
#     print('children', children)

    if (len(T.nodes()) == 1):
        combos[0] = {}

        return combos

    elif (sb):
        newNode = np.max(T.nodes()) + 1
        R = nx.Graph()
        R.add_edges_from(T.edges())
        R.add_edges_from([(newNode,h), (newNode,newNode+1), (newNode+1,h)])

        P = nx.Graph()
        P.add_edges_from(R.edges())

        aut = fnr_wc.findSubgraphInstances(R, P, False)
        tClasses = fnr_wc.findEquivalenceClasses(aut)

        reps = list(tClasses.keys())

#         print('reps:', reps)
#         print('neighbors of h', list(T[h]))

        children = []

        for node in reps:
            if (node in T[h]):
                children = children + [node] * len(tClasses[node])

#         print('children:', children)

        for i in range(len(children)):
            c = co.Counter()
            c.update(map(tuple, map(sorted, it.combinations(children, i+1))))

            l = list(c.keys())
            combo = []

            for x in l:
                x = np.array(x)
                y = np.array(x)

                for key in tClasses.keys():
                    if (key in T[h]):
                        eqList = list(tClasses[key])

                        indices = np.where(x==key)
#                         print('indices:', indices)

                        if (len(indices[0]) > 1):
                            replace = eqList[:len(indices[0])]
#                             print('replace:', replace)

                            for (ind, rep) in zip(indices[0], replace):
                                y[ind] = rep

                combo.append(tuple(y))

            combos[i+1] = combo

    else:
        for i in range(len(children)):
            combos[i+1] = list(it.combinations(children, i+1))

#     print('combos', h, ':', combos)
    return combos


# input: edges of a graph written in bracketed pairs
# output: graph with corresponding edges
def parseStringToGraph(s):
    l = []

    if(s[len(s)-1] == ' '):
        s = s[:-1]

    if(s[0] == '[' and s[len(s)-1] == ']'):
        s = s[-(len(s)-1):-1]
        start = [key for key, val in enumerate(s) if val == '[']
        end = [key for key, val in enumerate(s) if val == ']']
        comma = [key for key, val in enumerate(s) if (val == ',' and not s[key-1] == ']')]

        for i in range(len(start)):
            left = [v for k, v in enumerate(s) if k in range(start[i]+1, comma[i])]
            right = [v for k, v in enumerate(s) if k in range(comma[i]+1, end[i])]

            left = int(''.join(left))
            right = int(''.join(right))
            l.append((left, right))

    H = nx.Graph()
    H.add_edges_from(l)

    return H
