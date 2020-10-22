import numpy as np
import networkx as nx
import time
import pickle

import functions_regular as fnr
import functions_subtrees as fns

print('import done')

# this test script is going to compare the runtime of fnr and fns on a
# test network from the community_fitNet corpus. it will write the output
# to a file named test.out. the script will run on a single core of RMACC.

#########################################
########### helper functions ############
#########################################

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

def graphToString(G):
    edges = list(G.edges())
    string = '['

    for edge in edges:
        to_add = '[' + str(edge[0]) + ',' + str(edge[1]) + '], '
        string += to_add

    string = string[:len(string)-2]
    string += ']'

    return string

print('helper functions done')


#########################################
################ script #################
#########################################


networkFile = open('data/network48.txt', 'r')
networkString = networkFile.readlines()
networkString = networkString[0]
networkGraph = nx.Graph()
networkGraph = parseStringToGraph(networkString)

queryFile = open('data/query.txt', 'r')
queryString = queryFile.readlines()
queryString = queryString[0]
queryGraph = nx.Graph()
queryGraph = parseStringToGraph(queryString)

print('import data done')

start_fnr = time.time()
out_fnr = fnr.findSubgraphInstances(queryGraph, networkGraph, True)
delta_fnr = time.time() - start_fnr

print('regular done')

start_fns = time.time()
out_fns = fns.findSubgraphInstances(queryString, networkString)
delta_fns = time.time() - start_fns

print('subtrees done')

networkIndex = 48

if (out_fnr != out_fns):
    delta_fnr = -1
    delta_fns = -1
else:
    print('query:', queryString, 'network:', networkIndex, 'regular:', delta_fnr, 'subtrees:', delta_fns)
