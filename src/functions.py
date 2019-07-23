# this is for testing implementations of graph algorithms to speed up motif seach in a network.
# after refactoring, testing, verification, and new additions make this available as open-source.
# also, the comments have to be organized using a documentation tool. i like to throw comments
# everywhere becaue it gives consolation that i'm writing useful stuff.
# Tarek Tohme, June 2019

import numpy as np
import pandas as pd
import networkx as nx
import collections as co
import random as rnd
import permanent as per
import operator
import itertools as it

count = 0

##########################################################################
###################### useful functions and classes ######################
##########################################################################

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

    if(len(G[g]) >= len(H[h])):
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
def findEquivalenceClasses(aut):
    eq = {}

    for n in aut[0].getDomain():   # loop through all nodes in the domain (in order)
        E = set()

        for f in aut:
            E = E.union({f.applyMap(n)})   # gather all nodes that n can be mapped to in E

        if (E not in eq.values()):
            eq[n] = E

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


###################################################################################
###################### grochow-kellis motif search algorithm ######################
###################################################################################

# finds all instances of query graph H in network G
# H: query graph
# G: network to be queried
def findSubgraphInstances(H, G, withSBC=True):
    # print('##### FIND SUBGRAPH INSTANCES #####')
    J = nx.Graph()
    J.add_edges_from(G.edges())

    # instances of H found in G
    instances = []
    trees = 0

    if (len(H) > len(G)):
        return 0

    if(withSBC):

        K = nx.Graph()
        K.add_nodes_from(H.nodes())
        K.add_edges_from(H.edges())

        aut = findSubgraphInstances(H, K, False)  # returns list of automorphisms of H

        # print()
        # print('##################################################')
        # print('aut: ',end='')
        # print([i.getMap() for i in aut])

        M = symmetryConditions(aut)

        # print()
        #
        # print('M: ',end='')
        # print(M)

        eqClasses = findEquivalenceClasses(aut)
        HE = [t for t in eqClasses.keys()]

        # print('equivalence classes: ',end='')
        # print(eqClasses)

        # sort nodes of G by degree
        sortedDegreeG = sortDegrees(G, G.nodes())
        mapping = dict(zip(sortedDegreeG['node'], range(len(G))))
        G = nx.relabel_nodes(G, mapping)
        J = G

        onionDecompose(G)
        onionDecompose(H)

    else:
        HE = [i for i in H.nodes()]
        M = None


    for g in sorted(list(G.nodes())):
        for h in HE:

            # print('HE:', HE)
            # print('h:',h,'g:',g)
            # print('######### NEW NODE #########')

            if(canSupport(h, g, H, G)):
                r_g = False
                r_h = False

                if (withSBC):
                    r_g = G.node[g]['is_root']  #r_g is the child of g on g's tree, if g is a root
                    r_h = H.node[h]['is_root']

                if(r_g and r_h):
                    T = subtree(H, r_h, h)
                    T.add_edge(r_h, h)
                    S = subtree(G, r_g, g)
                    S.add_edge(r_g, g)

                    trees += countRootedSubtrees(T, h, S, g)

                else:
                    f = Map([(h, g)])  # initialize partial map with f(h) = g. argument must be list
                    iso = isomorphicExtensions(f, H, G, 1, M)  # M

                    # print([k.getMap() for k in iso])

                    if(type(iso) == type(f)):  # sometimes iso is single element
                        instances.append(iso)
                    else:
                        [instances.append(i) for i in iso]   # instances might contain duplicate maps

        G.remove_node(g)

    if(H.edges() == G.edges()):
        instances = bijectionsOnly(instances)

    if (withSBC):
        return len(instances) + trees

    return instances


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

    # check for induced isomorphism.
    for n in neighborsR:

        # print('partial map: ')
        # print(np.vstack(f.getMap()))
        #
        # print('try range extension (n): ',end='')
        # print(n)

        out = findCandidates(m, n, H, G, f)

        neighbMinD = out[0]
        f_neighbMinD = out[1]
        neighbNinR = out[2]

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
                # print('failed SBC')
                pass

        else:
            # print('failed isomorphism test')
            pass

    # for c in range(call):
    #     print('   ',end='')
    #
    # print('EXIT ISOMORPHIC EXTENSIONS CALL #', end='')
    # print(call)
    # for c in range(call):
    #     print('   ',end='')
    # print('output: ',end='')

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
# M[n] = S: S is a set of nodes such that l(n) < min(l(k) | k in S) i.e. l(S) must be > l(n)
def symmetryConditions(aut):
    M = {}  # dict containing nodes in H that have conditions and their constraining sets of nodes
    eqClasses = findEquivalenceClasses(aut)  # list of tuples: eq. class, representative node
    HE = [i for i in eqClasses.keys()]  # list of representative nodes

    # print(eqClasses)

    for n in HE:
        np = n
        A = aut
        temp_eqClasses = eqClasses

        while len(A) > 1:
            # print('A before:')
            # print([tuple(i.getMap()) for i in A])

            Sp = temp_eqClasses[np]
            Sp.remove(np)
            M[np] = Sp
            # print('M[',np,'] = ',Sp)

            A = [f for f in A if f.applyMap(np) == np]

            # print('A after:')
            # print([tuple(i.getMap()) for i in A])

            temp_eqClasses = findEquivalenceClasses(A)

            maxx = 0

            for k in temp_eqClasses.keys():
                size = len(temp_eqClasses[k])

                if (size > maxx):
                    maxx = size
                    np = k

    return M




###########################################################################
###################### onion decomposition algorithm ######################
###########################################################################

# labels each node with its onion layer and coreness, adds them as attributes to G
# G: network to be decomposed
def onionDecompose(G):
    K = nx.Graph()
    K.add_edges_from(G.edges())

    is_root = {}
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
            is_root[v] = None

            if (core == 2):
                shell_2.add(v)

            K.remove_node(v)  # delete from D
            count = count + 1

        layer = layer + 1

        D = [len(K[k]) for k in K.nodes()]

        if(D):
            minn = min(D)

            if (minn > core):
                core = minn

    for v in shell_2:
        for u in G[v]:
            if coreness[u] == 1:
                is_root[v] = u

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

    T = nx.Graph()
    T.add_node(0)

    for i in range(1, n):
        r = rnd.randint(0,i-1)
        T.add_edge(r, i)

    return T


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


# returns the subtree of T rooted at r with parent p
def subtree(T, r, p):
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
        p += per.permanent(A)

    return p


# draws a bipartite graph with labeled edges
def drawBipartiteGraph(F, l, r):

    M = nx.to_numpy_matrix(F)
    D = nx.to_dict_of_dicts(F)

    indices = [i for i in range(l)]
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
def countMatchings(T, rootT, S, rootS, call=0, init=True):

    parentT = rootT
    parentS = rootS
    childrenT = set(T[rootT]).difference({rootT})
    childrenS = set(S[rootS]).difference({rootS})
    G = nx.Graph()

    # base case
    if (not childrenT):
        return 1

    if (len(childrenT) > len(childrenS)):
        return 0

    indices = [0]

    for key_i, i in enumerate(childrenT):
        for key_j, j in enumerate(childrenS):
            Ti = subtree(T, i, rootT)
            Sj = subtree(S, j, rootS)
            call2 = call + 1
            k = countMatchings(Ti, i, Sj, j, call2, False)

            if (k > 0):
                G.add_edge(key_i, key_j+len(childrenT), weight = k)
        indices.append(key_i)

    l = len(childrenT)
    r = len(childrenS)

    M = nx.to_numpy_matrix(G, [i for i in range(l+r)])

    matchings = rectPermanent(M, l, r)

    return matchings


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

    # recursion
    count = countMatchings(T, rootT, S, rootS)
    return count
