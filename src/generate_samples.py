# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


# To run, type:
# python3 generate_samples.py graph_file_name sample_file_name sample_size data_directory [random_seed]
# For example:
# python3 generate_samples.py graph_random_bowfree_10nodes-0.txt sample-100.txt 100 ../data/test/


import numpy as np
import networkx as nx
from typing import Sequence, Iterable
import sys
from read_graph_from_file import read_graph_from_file


def simulate_sem_multivariate_gaussian(
        D: Iterable[Sequence[int]],
        B: Iterable[Sequence[int]],
        n: int) -> np.ndarray:
    """Simulate samples from SEM with specified type of noise.

    Args:
        D: adjacency matrix of directed edges
        B: adjacency matrix of bidirected edges
        n: number of samples

    Returns:
        X: [n,d] sample matrix
        delta: [d,d]
        beta: [d,d]
    """
    G = nx.DiGraph(np.array(D))
    W = nx.to_numpy_array(G)
    d = W.shape[0]
    x_dims = 1
    X = np.zeros([n, d, x_dims])
    ordered_vertices = list(nx.topological_sort(G))
    assert len(ordered_vertices) == d
    delta = np.zeros([d, d])
    beta = np.zeros([d, d])

    low = 0.5
    high = 2.0

    diff = high-low
    for i in range(d):
        for j in range(d):
            if D[i][j] > 0:
                delta[i][j] = np.random.uniform(-diff, diff)
                if delta[i][j] < 0.0:
                    delta[i][j] -= low
                else:
                    delta[i][j] += low

    low = 1.0
    high = 2.0

    diff = high-low
    for i in range(d-1):
        for j in range(i+1, d):
            if B[i][j] > 0:
                beta[i][j] = np.random.uniform(-diff, diff)
                if beta[i][j] < 0.0:
                    beta[i][j] -= low
                else:
                    beta[i][j] += low
                beta[j][i] = beta[i][j]
    
    low = 0.7
    high = 1.2

    for i in range(d):
        sum = 0.0
        for j in range(d):
            if i != j:
                sum += abs(beta[i][j])

        beta[i][i] = np.random.uniform(low, high)
        if beta[i][i] < 0.0:
            beta[i][i] -= (low + sum)
        else:
            beta[i][i] += (low + sum)

    epsilon = np.random.multivariate_normal([0]*d,beta,n)
    for j in ordered_vertices:
        parents = list(G.predecessors(j))
        eta = X[:, parents, 0].dot(delta[parents, j])
        X[:, j, 0] = eta + epsilon[:,j]

    return X, delta, beta


def main():

    data_directory = ''
    graph_file_name = ''
    sample_file_name = ''
    if len(sys.argv) > 1:
        graph_file_name = sys.argv[1]
    if len(sys.argv) > 2:
        sample_file_name = sys.argv[2]
    sample_size = 1000
    if len(sys.argv) > 3:
        sample_size = int(sys.argv[3])
    if len(sys.argv) > 4:
        data_directory = sys.argv[4]
    random_seed = 12345
    if len(sys.argv) > 5:
        random_seed = int(sys.argv[5])
    np.random.seed(random_seed)
    delta_beta_file_name = data_directory + 'delta_beta_' + graph_file_name
    graph_file_name = data_directory + graph_file_name
    sample_file_name = data_directory + sample_file_name
        
    graph = read_graph_from_file(graph_file_name)
    data, delta, beta = simulate_sem_multivariate_gaussian(graph[0], graph[1], sample_size)

    data = data[:,:,0]

    nrows,ncols = data.shape
    wrt = ''
    for i in range(nrows):
        for j in range(ncols-1):
            wrt = wrt + str(data[i][j]) + ' '
        wrt = wrt + str(data[i][ncols-1]) + '\n'
    f = open(sample_file_name,"w")
    f.write(wrt)
    f.close()

    wrt = 'nnodes\n'
    wrt = wrt + str(ncols) + '\n'
    wrt = wrt + 'delta\n'
    for i in range(ncols):
        for j in range(ncols-1):
            wrt = wrt + str(delta[i][j]) + ' '
        wrt = wrt + str(delta[i][ncols-1]) + '\n'
    wrt = wrt + 'beta\n'
    for i in range(ncols):
        for j in range(ncols-1):
            wrt = wrt + str(beta[i][j]) + ' '
        wrt = wrt + str(beta[i][ncols-1]) + '\n'
    f = open(delta_beta_file_name,"w")
    f.write(wrt)
    f.close()

    

if __name__ == "__main__":

    main()
