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


# returns only functions that are unique in the set maps
# maps: set of maps, possibly with duplicates
def returnUniqueMaps(maps):
    dummy = []
    to_return = []

    for f in maps:
        if (f.isBijection()):
            dummy.append(f)

    for g in dummy:
        if(g.inSet(to_return)):
            pass

        else:
            to_return.append(g)

    return to_return


###############################################
#### grochow-kellis motif search algorithm ####
###############################################

# finds all instances of query graph H in network G
# H: query graph
# G: network to be queried
def findSubgraphInstances(H, G, withSBC=True):
    # print('##### FIND SUBGRAPH INSTANCES #####')

    # instances of H found in G
    instances = []

    if(withSBC):
        aut = findSubgraphInstances(H, H, False)  # returns list of automorphisms of H
        M = symmetryConditions(aut)
        print(M)

        eqClasses = findEquivalenceClasses(aut)
        HE = [t[1] for t in eqClasses]

    else:
        HE = H
        M = None

    # sort nodes of G by degree
    sortedDegreeG = sortDegrees(G, G.nodes())

    for g in sortedDegreeG['node']:
        # print(g)
        for h in HE:
            # print('######### NEW NODE #########')
            if(canSupport(h, g, H, G)):
                f = Map([(h, g)])  # initialize partial map with f(h) = g. argument must be list
                iso = isomorphicExtensions(f, H, G, 1, M)  # M
                # print([k.getMap() for k in iso])

                if(iso):  # sometimes iso is empty
                    [instances.append(i) for i in iso]   # instances might contain duplicate maps

        # G.remove_node(g)

    # instances = [next(iter(k)) for k in instances]  # removes the elements from their sets

    to_return = returnUniqueMaps(instances)

    return to_return


# finds all isomorphic extensions of a partial map [satisfying the symmetry-breaking
# condition C at h]. returns them in a list of tuples [(a1, b1), ... ,(ak, bk)]
# f: partial map to be extended
def isomorphicExtensions(f, H, G, call, M = None): # M
    # for c in range(call):
    #     print('   ',end='')
    #
    # print('ISOMORPHIC EXTENSIONS CALL #',end='')
    # print(call)


    isomorphisms = []
    neighborsR = {None}
    neighborsD = {None}
    D = f.getDomain()
    R = f.getRange()

    if(set(D) == set(H.nodes())):
        # for c in range(call):
        #     print('   ',end='')
        # print('INSTANCE FOUND ON CALL #', end='')
        # print(call)

        f.extend(list(np.sort(f.getMap(), order='domainNode')))
        # for c in range(call):
        #     print('   ',end='')
        # print(list(f.getMap()))
        return f

    m = mostConstrainedNode(D, H)
    # print('domain extension (m): ', end='')
    # print(m)

    # list of neighbors of f(D)
    for r in R:
        neighborsR = neighborsR.union(set(G[r]))

    # list of neighbors of D
    for d in D:
        neighborsD = neighborsD.union(set(H[d]))

    # exclusive neighborhood
    neighborsR = neighborsR.difference(set(R)).difference({None})
    neighborsD = neighborsD.difference(set(D)).difference({None})

    # print('neighbors of partial range: ',end='')
    # print(neighborsR)
    # print('neighbors of partial domain: ',end='')
    # print(neighborsD)

    # check for induced isomorphism
    for n in neighborsR:
        # print(chr(n+65), end='')

        # print('partial map: ')
        # print(np.vstack(f.getMap()))
        #
        # print('try range extension (n): ',end='')
        # print(n)

        neighbMinD = set(H[m]).intersection(set(f.getDomain()))
        f_neighbMinD = set([f.applyMap(k) for k in neighbMinD])
        neighbNinR = set(G[n]).intersection(set(f.getRange()))

        # print('neighbors of m in D: ', end='')
        # print(neighbMinD)
        # print('neighbors of n in R: ', end='')
        # print(neighbNinR)
        # print('f(neighbors of m in D): ', end='')
        # print(f_neighbMinD)

        if(set(f_neighbMinD) == set(neighbNinR)): # add condition for tree string matching later

            if(checkSBC(m, n, M, f) or (M == None)):  # n conforms to symmetry-breaking conditions
                # print('range extension valid')
                fp = Map(list(f.getMap()))
                newNode = np.array((m, n), dtype=[('domainNode', 'i4'), ('rangeNode', 'i4')])

                fp.extend([newNode])
                # print('called')
                call2 = call + 1
                iso = isomorphicExtensions(fp, H, G, call2, M)
                # print('returned')

                if(type(iso) == type(fp)):
                    isomorphisms.append(iso)
                else:
                    [isomorphisms.append(i) for i in iso]

            else:
                pass

        else:
            # print('range extension valid')
            pass

    # for c in range(call):
    #     print('   ',end='')
    #
    # print('EXIT ISOMORPHIC EXTENSIONS CALL #', end='')
    # print(call)
    # for c in range(call):
    #     print('   ',end='')
    # print('output: ',end='')
    #
    # print([tuple(k.getMap()) for k in isomorphisms])

    return isomorphisms


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

    for key, value in M.items():
        if (m == key):

            if (value):  # case 1
                images = [f.applyMap(k) for k in value] # {L(k)|k in values}

                if (n < min(images)): # all values mapped, n conforms to sbc
                    # print('conforms to sbc. node accepted')
                    return True

                elif(min(images) == -1): # some values ont mapped
                    # print('some conditions unknown. node accepted')
                    return True

                else:  # n violates sbc
                    # print('violates sbc. node rejected')
                    return False

            else:  # case 4
                # print('no conditions on n. node accepted')
                return True

        elif (m in value):   # case 2

            if (f.applyMap(key) > 0):

                if (f.applyMap(key) < n):
                    # print('conforms to sbc. node accepted')
                    return True

                else:
                    # print('violates sbc. node rejected')
                    return False

            else:
                # print('some conditions unknown. node accepted')
                return True


    inValues = [m in k for k in M.values()]

    if (sum(inValues) == 0 and (m not in M.keys())):
        # print('no conditions on n. node accepted')
        return True

    print('something went wrong')
    return False


#######################################
#### onion decomposition algorithm ####
#######################################

# labels each node with its onion layer and coreness
# G: network to be decomposed
def onionDecompose(G):
    labels = {}
    core = 1
    layer = 1

    while (len(G.nodes()) > 0):
        thisLayer = [v for v in G.nodes() if len(G[v]) <= core]

        for v in thisLayer:
            newLabel = {'coreness': core,
                        'layer': layer}

            labels[v] = newLabel
            G.remove_node(v)  # delete from D

        layer = layer + 1

        D = [len(G[k]) for k in G.nodes()]

        if(D):
            minn = min(D)

            if (minn > core):
                core = minn

    return labels
