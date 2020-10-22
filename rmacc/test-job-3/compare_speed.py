import numpy as np
import networkx as nx
import time
import pickle
import sys

import functions_regular_wc as fnr_wc
import functions_subtrees_wc as fns_wc

# print('import done')

# this script runs all the speed tests we'll need for the paper on a few
# sample networks and query graphs, using the load balancer utility

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

# print('helper functions done')

#########################################
################ script #################
#########################################

# take the network index as a command line argument and run the querys on that network

networkIndex = int(sys.argv[1])
querySize = sys.argv[2]
queryStartIndex = int(sys.argv[3])      #inclusive
queryEndIndex = int(sys.argv[4])        #exclusive

f = open('data/CommunityFitNet_updated.pickle', 'rb')
data = pickle.load(f)
data_edgelists = data['edges_id']

g = open('data/allgraphs'+querySize+'indexed.txt', 'r')
query_data = g.readlines()

networkGraph = nx.Graph()
edges = data_edgelists.iloc[networkIndex]
networkGraph.add_edges_from(edges)
networkString = graphToString(networkGraph)

networkSize = len(networkGraph)

output = {}

# all query graphs of size k
for key, val in enumerate(query_data[queryStartIndex:queryEndIndex]):

    spaceIndex = val.find(' ')
    queryKey = val[:spaceIndex]
    queryString = val[spaceIndex+1:].rstrip()

    start_time_fnr = time.time()
    out_fnr = fnr_wc.findSubgraphInstances(queryString, networkString, True)
    delta_fnr = time.time() - start_time_fnr

    start_time_fns = time.time()
    out_fns = fns_wc.findSubgraphInstances(queryString, networkString)
    delta_fns = time.time() - start_time_fns

    if (out_fnr != out_fns):
        output[networkIndex] = {queryKey: 'error'}
    else:
        output[networkIndex] = {queryKey: {'runtime_fnr': delta_fnr,
                                          'runtime_fns': delta_fns}}

print(output)
