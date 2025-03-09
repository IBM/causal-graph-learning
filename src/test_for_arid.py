# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


# To run, type:
# python3 test_for_arid.py


import itertools


def depth_first_search_traversal(diedges, node, visited):
    visited.add(node)

    for j in range(len(diedges[node])):
        head = diedges[node][j]
        if head not in visited:
            depth_first_search_traversal(diedges, head, visited)


def does_subgraph_contain_ctree(graph, dist_subset):
    D = graph['di']
    B = graph['bi']
    nnodes = len(D)

    odiedges = [[] for i in range(nnodes)]
    idiedges = [[] for i in range(nnodes)]
    biedges = [[] for i in range(nnodes)]

    for i in range(nnodes):
        if i in dist_subset:
            for j in range(nnodes):
                if j in dist_subset:
                    if D[i][j] == 1:
                        odiedges[i].append(j)
                        idiedges[j].append(i)
                    if B[i][j] == 1:
                        biedges[i].append(j)

    # check that every node has incident bidirected edge
    for node in dist_subset:
        if len(biedges[node]) == 0:
            return False

    # check that all but one nodes have positive outdegree
    root = None
    for node in dist_subset:
        if len(odiedges[node]) == 0:
            if root == None:
                root = node
            else:
                return False # there can't be more than one root node

    # check that root node exists (root must have positive indegree)
    if root == None or len(idiedges[root]) == 0:
        return False

    # check that arborescence exists
    visited = set()
    depth_first_search_traversal(idiedges, root, visited)
    for node in dist_subset:
        if node not in visited:
            return False

    return True

def does_graph_contain_ctree(graph):
    D = graph['di']
    B = graph['bi']
    nnodes = len(D)

    distnodes = []
    for i in range(nnodes):
        for j in range(nnodes):
            if B[i][j] == 1:
                distnodes.append(i)
                break

    ndistnodes = len(distnodes)
    for k in range(4,ndistnodes+1):
        # find all subsets of distnodes of size k
        all_subsets = [list(i) for i in itertools.combinations(distnodes,k)]
        # iterate over all subsets, add relevant parents and see if arborescence exists
        for i in range(len(all_subsets)):
            subgraph_contains_ctree = does_subgraph_contain_ctree(graph, all_subsets[i])
            if subgraph_contains_ctree == True:
                return True

    return False


def is_there_a_cycle(diedges, node, visited, inpath):
    visited[node] = True
    inpath[node] = True

    for j in range(len(diedges[node])):
        head = diedges[node][j]
        if visited[head] == False:
            graph_has_cycle = is_there_a_cycle(diedges, head, visited, inpath)
            if graph_has_cycle == True:
                return True
        elif inpath[head] == True:
            return True

    inpath[node] = False
    return False

def is_directed_graph_acyclic(D):
    nnodes = len(D)
    diedges = [[] for i in range(nnodes)]
    for i in range(nnodes):
        for j in range(nnodes):
            if D[i][j] == 1:
                diedges[i].append(j)
    
    visited = [False for i in range(nnodes)]
    inpath = [False for i in range(nnodes)]
    for node in range(nnodes):
        if visited[node] == False:
            graph_has_cycle = is_there_a_cycle(diedges, node, visited, inpath)
            if graph_has_cycle == True:
                return False
    return True

def is_graph_bowfree(graph):
    D = graph['di']
    B = graph['bi']
    nnodes = len(D)

    for i in range(nnodes-1):
        for j in range(i+1, nnodes):
            if B[i][j] == 1 and (D[i][j] == 1 or D[j][i] == 1):
                return False

    return True


def is_graph_arid(graph):

    D = graph['di']
    B = graph['bi']
    nnodes = len(D)
    
    graph_is_bowfree = is_graph_bowfree(graph)
    if graph_is_bowfree == False:
        return False

    graph_is_acyclic = is_directed_graph_acyclic(D)
    if graph_is_acyclic == False:
        return False

    nnodesdist = 0 # number of nodes in district
    for i in range(nnodes):
        for j in range(nnodes):
            if B[i][j] == 1:
                nnodesdist += 1
                break

    # only c-components with 4 or more nodes in the district can be non-arid
    if nnodesdist < 4:
        return True

    graph_contains_ctree = does_graph_contain_ctree(graph)
    if graph_contains_ctree == True:
        return False
    
    return True



def test_for_arid(nodes, parents, edges):
    # nodes contains the nodes in the district
    # there might be nodes in parents that are not in the district
    all_nodes = []
    for node in nodes:
        all_nodes.append(node)
    for i in range(len(nodes)):
        for j in range(len(parents[i])):
            node = parents[i][j]
            if node not in all_nodes:
                all_nodes.append(node)
    all_nodes.sort()

    nnodes = len(all_nodes)
    D = [[0]*nnodes for i in range(nnodes)]
    B = [[0]*nnodes for i in range(nnodes)]

    # since nodes might not be numbered from 0 to nnodes-1,
    # we are going to map the real node numbers to
    # the range 0 to nnodes-1
    for i in range(len(nodes)):
        for j in range(len(parents[i])):
            head = nodes[i]
            index_head = all_nodes.index(head)
            tail = parents[i][j]
            index_tail = all_nodes.index(tail)
            D[index_tail][index_head] = 1

    for i in range(len(edges)):
        index1 = all_nodes.index(edges[i][0])
        index2 = all_nodes.index(edges[i][1])
        B[index1][index2] = 1
        B[index2][index1] = 1

    graph = {}
    graph['di'] = D
    graph['bi'] = B

    graph_is_arid = is_graph_arid(graph)
    return graph_is_arid
        

def main():
    nodes = (0, 1, 2, 3, 5, 6, 8, 9)
    edges = ((0, 1), (0, 3), (1, 2), (1, 3), (1, 6), (1, 9), (2, 8), (3, 6), (5, 6))
    parents = ((), (7,), (0,), (2, 7), (0, 2, 4), (4,), (0, 4, 6, 7), (0, 2, 3, 4, 5, 7))
    if test_for_arid(nodes,parents,edges) == True:
        print('c-comp with nodes '+str(nodes)+', edges '+str(edges)+', parents '+str(parents)+', is arid')
    else:
        print('c-comp with nodes '+str(nodes)+', edges '+str(edges)+', parents '+str(parents)+', is not arid')


if __name__ == "__main__":

    main()
