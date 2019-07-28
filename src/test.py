import numpy as np
import pandas as pd
import networkx as nx
import functions_regular as fnr
import functions_subtrees as fns

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


# reads the text file output by Josh and returns test name, hstring,
# gstring and answer in a dict
def parseFileToDict(file):
    f = open(file, 'r')
    D = {}
    index = 0
    first = True

    for line in f:
        if (line[0:4] == 'Test'):
            title = line[:-1]

        if (line[0] == '['):
            if (first):
                hstring = line[:-1]
                first = False
            else:
                gstring = line[:-1]
                first = True


        if (ord(line[0]) >= 48 and ord(line[0]) <= 57):
            ans = int(line)
            D[index] = {'title': title,
                        'hstring': hstring,
                        'gstring': gstring,
                        'answer': ans}
            index = index + 1

    return D


def test_findSubgraphInstances(file, which='regular'):

    D = parseFileToDict(file)

    for key in D.keys():
        print(D[key]['title'])
        H = nx.Graph()
        G = nx.Graph()

        hGraph = parseStringToGraph(D[key]['hstring'])
        gGraph = parseStringToGraph(D[key]['gstring'])

        H.add_edges_from(hGraph.edges())
        G.add_edges_from(gGraph.edges())

        if (which == 'regular'):
            count = fnr.findSubgraphInstances(H, G, True)

        if (which == 'subtrees'):
            count = fns.findSubgraphInstances(H, G, True)

        l = count
        ans = D[key]['answer']

        print('H: ',end='')
        print(D[key]['hstring'])
        print('G: ',end='')
        print(D[key]['gstring'])
#         print('instances: ',end='')
#         print([tuple(i.getMap()) for i in instances])
        print('correct number of instances: ',end='')
        print(ans)
        print('instances found: ',end='')
        print(l)

        if(l == ans):
            print('succeeded')
        else:
            print('failed')

        print()
