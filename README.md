# Integer Programming Based Methods and Heuristics for Causal Graph Learning

This repository contains the code and the data used for the experiments in the paper
"Integer Programming Based Methods and Heuristics for Causal Graph Learning" by Sanjeeb Dash, Joao Goncalves, Tian Gao. The paper was accepted at AISTATS 2025.


## Requirements:
* The code was tested only on Linux.
* The code is written in python and was tested with python 3.9.21.
* The code uses the following python packages: numpy, networkx, scipy. It was tested with numpy 1.26.4, networkx 3.2.1, scipy 1.12.0.
* The code uses the commercial solver [IBM ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio). It was tested with version 22.1.0.0.


## Folders description:
* `src`: contains the code.
* `data`: contains the datasets used in the paper.
* `results`: contains the output graphs from the experiments in the paper.


## Installation:
* Install the following python packages: numpy, networkx, scipy.
* Install IBM ILOG CPLEX and setup the python interface.


## How to run code to learn graph:
`python3 learn.py -s score_file_name.pkl -i path_input_directory/ -o path_output_directory/ -h heuristics -c cuts -t MIP_time_limit_(seconds) -g type_of_graph_desired > log_file.txt`\
where:\
`score_file_name.pkl` = file with input c-components and corresponding scores,\
`path_input_directory/` = path to the directory where `score_file_name.pkl` resides,\
`path_output_directory/` = path to the directory where the output is written,\
`heuristics` = string that can be any combination of the letters a, b, c (i.e., a, b, c, ab, ac, bc, abc); a corresponds to operation (1) in the paper, b corresponds to operation (2) in the paper, and c corresponds to operation (3) in the paper,\
`cuts` = string that can be any combination of the letters `a`, `b`, `c` (i.e., `a, b, c, ab, ac, bc, abc`); `a` corresponds to theorem 3.2 in the paper applied to a 3-node graph and the set S being a pair of nodes, `b` corresponds to theorem 3.2 in the paper applied to a 3-node graph and the set S being the 3 nodes, and `c` corresponds to theorem 3.3 in the paper applied to a 3-node graph and the set S being the 3 nodes,\
`MIP_time_limit_(seconds)` = integer with the time limit for the MIP,\
`type_of_graph_desired` = either `aadmg`, `arid`, or `bowfree`.\
\

Example of a run with bowfree:\
`python3 learn.py -s score_sample_graph_random_bowfree_10nodes-0-cs2_sps3_ocps2.pkl -i ../data/test/ -o ../results/test/ -t 7200 -g bowfree > log_file.txt`\
\

Example of a run with bowfree and cuts from theorem 3.2:\
`python3 learn.py -s score_sample_graph_random_bowfree_10nodes-0-cs2_sps3_ocps2.pkl -i ../data/test/ -o ../results/test/ -cuts ab -t 7200 -g bowfree > log_file.txt`\
\

Example of a run with bowfree and heuristics `a`, `b`, and `c`:\
`python3 learn.py -s score_sample_graph_random_bowfree_10nodes-0-cs1_sps3_ocps1.pkl -i ../data/test/ -o ../results/test/ -h abc -t 7200 -g bowfree > log_file.txt`\


## How to generate random bow-free graphs:
`python3 generate_random_graph.py ngraphs nnodes pdir pbidir graph_file_name [graph_directory/]`\
where:\
`ngraphs` = number of random graphs to generate\
`nnodes` = number of nodes in each graph\
`pdir` = probability of generating a directed edge\
`pbidir` = probability of generating a bidirected edge\
`graph_file_name` = base name for the graphs\
`graph_directory` = (optional) directory where the graph files will be written\
\
Example:\
`python3 generate_random_graph.py 2 5 0.4 0.3 graph_random ../data/test/`\


## How to generate samples:
`python3 generate_samples.py graph_file_name sample_file_name sample_size data_directory [random_seed]`\
where:\
`graph_file_name` = name of file containing graph\
`sample_file_name` = name of file where samples are going to be written\
`sample_size` = the number of samples to generate\
`data_directory` = the directory where the graph_file_name exists and where the file sample_file_name will be written\
`random_seed` = (optional) random seed\
\
Example:\
`python3 generate_samples.py graph_random_bowfree_10nodes-0.txt sample-100.txt 100 ../data/test/`\


## How to generate initial c-components and corresponding scores:
`python3 generate_ccomponents.py sample_file_name max_district_size max_parent_set_size_1 max_parent_set_size_2 data_directory`\
where:\
`sample_file_name` = name of file where samples are going to be written\
`max_district_size` = maximum number of nodes in a district\
`max_parent_set_size_1` = maximum number of nodes in the parent set of the node in the one-node district\
`max_parent_set_size_2` = maximum number of nodes in the parent set of each node in the district (district with 2 or more nodes)\
`data_directory` = the directory where the sample_file_name exists and where the pickle file with the c-components and scores will be written\
The output file will have the following name: `score_sample_file_name-csX-spsY-ocpsZ.pkl`, where `X` is max_district_size, `Y` is max_parent_set_size_1, and `Z` is max_parent_set_size_2.\
\
Example:\
`python3 generate_ccomponents.py sample_graph_random_bowfree_10nodes-1.txt 2 3 2 ../data/test/`\

## How to compute the total score of a graph and the scores of the individual c-components:
`python3 evaluate_scores.py graph_file_name sample_file_name data_directory`\
where:\
`graph_file_name` = name of file containing graph\
`sample_file_name` = name of file containing samples\
`data_directory` = the directory where the graph_file_name and sample_file_name exist\

Example:\
`python3 evaluate_scores.py graph_random_bowfree_10nodes-0.txt sample_graph_random_bowfree_10nodes-0.txt ../data/test/`\
