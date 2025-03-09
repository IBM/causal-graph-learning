# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


# To run, type:
# python3 generate_random_graph.py ngraphs nnodes pdir pbidir graph_file_name [data_directory]
# The files names are graph_file_name-i.txt, where i is the number of the graph
# For example:
# python3 generate_random_graph.py 2 5 0.4 0.3 graph_random ../data/test/


import random
import sys

def write_graph_to_file(graph, filename):
    """ Writes a graph to a file

    The following is an example of an output file:

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

    file = open(filename, 'w')
    file.write('n_nodes n_directed n_bidirected\n')

    D = graph[0]
    B = graph[1]
    nnodes = len(D)

    ndir = 0
    for i in range(nnodes):
        for j in range(nnodes):
            if D[i][j] == 1:
                ndir += 1

    nbidir = 0
    for i in range(nnodes-1):
        for j in range(i+1, nnodes):
            if B[i][j] == 1:
                nbidir += 1

    file.write(str(nnodes) + ' ' + str(ndir) + ' ' + str(nbidir) + '\n')

    file.write('directed\n')
    for i in range(nnodes):
        for j in range(nnodes):
            if D[i][j] == 1:
                file.write(str(i) + ' ' + str(j) + '\n')

    file.write('bidirected\n')
    for i in range(nnodes-1):
        for j in range(i+1, nnodes):
            if B[i][j] == 1:
                file.write(str(i) + ' ' + str(j) + '\n')
    
    file.close()

    
def generate_random_bowfree_graph(nnodes, pdir, pbidir):

    D = [[0]*nnodes for i in range(nnodes)]
    B = [[0]*nnodes for i in range(nnodes)]
    
    for i in range(nnodes-1):
        for j in range(i+1, nnodes):
            if random.random() < pdir:
                D[i][j] = 1
            elif random.random() < pbidir:
                B[i][j] = 1
                B[j][i] = 1
    
    graph = {}
    graph[0] = D
    graph[1] = B
    return graph

def remove_bidirected_edges(graph):
    # each node can touch at most one bidirected edge
    D = graph[0]
    B = graph[1]
    nnodes = len(D)
    
    for j in range(nnodes):
        dipar = []
        bipar = []
        for i in range(nnodes):
            if B[i][j] == 1:
                bipar.append(i)

        n_in_edges = len(bipar)
        while n_in_edges > 1:
            index = random.randint(0, n_in_edges-1)
            B[bipar[index]][j] = 0
            B[j][bipar[index]] = 0
            bipar.pop(index)
            n_in_edges -= 1
    
    return graph

def remove_edges_until_max_in_arrows(graph, max_in_arrows):
    # max_in_arrows is the maximum number of arrows pointing to each node
    D = graph[0]
    B = graph[1]
    nnodes = len(D)
    
    for j in range(nnodes):
        dipar = []
        bipar = []
        for i in range(nnodes):
            if D[i][j] == 1:
                dipar.append(i)
            if B[i][j] == 1:
                bipar.append(i)

        n_in_edges = len(dipar)+len(bipar)
        while n_in_edges > max_in_arrows:
            index = random.randint(0, n_in_edges-1)
            if index < len(dipar):
                D[dipar[index]][j] = 0
                dipar.pop(index)
            else:
                index -= len(dipar)
                B[bipar[index]][j] = 0
                B[j][bipar[index]] = 0
                bipar.pop(index)
            n_in_edges -= 1
    
    return graph

    
def main():
    random.seed(1)
    data_directory = ''
    ngraphs = None # number of graphs to generate
    nnodes = None # number of nodes of the graph we want to generate
    pdir = None # probability of directed edge
    pbidir = None # probability of bidirected edge
    graph_file_name = None # filename where graph is going to be written
    if len(sys.argv) > 1:
        ngraphs = int(sys.argv[1])
    if len(sys.argv) > 2:
        nnodes = int(sys.argv[2])
    if len(sys.argv) > 3:
        pdir = float(sys.argv[3])
    if len(sys.argv) > 4:
        pbidir = float(sys.argv[4])
    if len(sys.argv) > 5:
        graph_file_name = sys.argv[5]
    if len(sys.argv) > 6:
        data_directory = sys.argv[6]
    graph_file_name = data_directory + graph_file_name

    for i in range(ngraphs):
        graph = generate_random_bowfree_graph(nnodes, pdir, pbidir)
#        graph = remove_bidirected_edges(graph)
#        max_in_arrows = 3
#        graph = remove_edges_until_max_in_arrows(graph, max_in_arrows)
        write_graph_to_file(graph, graph_file_name + '-' + str(i) + '.txt')


if __name__ == "__main__":

    main()
