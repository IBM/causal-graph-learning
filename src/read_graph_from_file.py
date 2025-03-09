# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


def read_graph_from_file(filename):
    """ Reads a graph from a file

    The following is an example of an input file:

    n_nodes n_directed n_bidirected
    4 3 2
    directed
    0 2
    2 3
    3 1
    bidirected
    0 1
    0 3

    """

    file = open(filename, 'r')
    data = file.readlines()
    file.close()

    # line 0 contains the heading: n_nodes n_directed n_bidirected
    line = data[1] # line 1 contains the values of n_nodes n_directed n_bidirected
    sline = line.strip('\n').split(' ')
    nnodes = int(sline[0])
    ndir = int(sline[1])
    nbidir = int(sline[2])

    D = [[0]*nnodes for i in range(nnodes)]
    B = [[0]*nnodes for i in range(nnodes)]

    # line 2 contains the heading: directed
    # the following lines contain the directed edges
    for nline in range(3,3+ndir):
        line = data[nline]
        sline = line.strip('\n').split(' ')
        tail = int(sline[0])
        head = int(sline[1])
        D[tail][head] = 1

    # line 3+ndir contains the heading: bidirected
    # the following lines contain the bidirected edges
    start = 3+ndir+1
    for nline in range(start,start+nbidir):
        line = data[nline]
        sline = line.strip('\n').split(' ')
        tail = int(sline[0])
        head = int(sline[1])
        B[tail][head] = 1
        B[head][tail] = 1

    graph = {}
    graph[0] = D
    graph[1] = B
    return graph
