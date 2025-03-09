# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


from scipy.special import comb
from utils import score_write_to_file, compute_BiC_c_comps_sets, check_connected_component,\
    find_bi_connected_node, has_cycle_in_c_comp
import itertools


def generate_scores_bidirect(data, single_c_parent_size, other_c_parent_size, c_size, file_name):
    # generate data for all bidirected scores
    num_sample, num_var = data.shape

    # compute scores for each graph
    scores = {}
    all_c_components = gererate_all_c_components(num_var, c_size)

    for iComp_size in all_c_components:
        c_comps_of_size = all_c_components[iComp_size]
        parent_size = single_c_parent_size if iComp_size < 2 else other_c_parent_size
        # ------------------
        # this is to generate c-comps with more than 2 node districts without parents
        #if iComp_size > 2:
        #    parent_size = 0
        # ------------------
        for each_comp in c_comps_of_size: # in a list

            edges = gererate_all_bidirected_edges(each_comp)
            for each_edge_set in edges:

                if iComp_size > 1:
                    save_edge_set = each_edge_set
                else:
                    save_edge_set = ()

                dic_key = (each_comp, save_edge_set)
                scores[dic_key] = {}

                parent_set_list = []
                for each_var_in_comp in each_comp:
                    # generate parental size
                    parent_set_list.append( gererate_all_parent_sets(each_var_in_comp,
                                                                     each_comp,
                                                                     num_var,
                                                                     parent_size,
                                                                     save_edge_set) )

                # generate parent set combination
                all_parent_list = list(itertools.product(*parent_set_list))

                # compute scores
                for each_parent_set_config in all_parent_list:

                    if not has_cycle_in_c_comp(each_parent_set_config, each_comp):

                        scores[dic_key][each_parent_set_config] = compute_BiC_c_comps_sets(data,
                                                                                            each_comp,
                                                                                            each_parent_set_config,
                                                                                           each_edge_set)

    # save files
    print('saving file..\n')
    score_write_to_file([data, scores], file_name)

    return


def gererate_all_bidirected_edges(each_comp):
    # generate all possible bi-directed edges within a c component
    edge_num = len(each_comp) - 1
    bi_edges = []
    edge_num_upper_bound = comb(len(each_comp), 2, exact=True)

    if edge_num < 1:
        flatten_list = [-1]  # no bidrected edges

    else:
        all_edges = list(itertools.combinations( sorted(set(each_comp)), 2))
        # for iCsize in range(1, edge_num+1): # size 0 to c_size
        for iCsize in range(edge_num, edge_num_upper_bound+1):  # size edge_num to max

            sets = list(itertools.combinations( set(all_edges), iCsize))

            # check connected components, make sure nodes are connected to each other in c-comp
            sets = check_connected_component(sets, each_comp)
            bi_edges.append(sets)

        flatten_list = [item for sublist in bi_edges for item in sublist]
    return flatten_list

def gererate_all_parent_sets(iVar,
                             var_in_comp,
                             num_var,
                             c_size,
                             comp_edges):
    # generate all possible c components
    # comp_edges: edges in a comp

    # bi_directed_edge: bi directed edges from iVar
    bi_directed_node = find_bi_connected_node(iVar, comp_edges)
    candidate_parent_in_comp =  set(list(var_in_comp)) - set([iVar])- \
                                set(list(bi_directed_node))
    c_comp = []
    for iCsize in range(c_size+1): # size 0 to c_size

        # add candiate parental sets
        parent_candidates =  list(sorted(
            (set(list(range(num_var))) - set(var_in_comp)).union(candidate_parent_in_comp)))
        # parent_candidates.remove(iVar)
        sets = list(itertools.combinations( sorted(set(parent_candidates)), iCsize))
        c_comp.append(sets)

    flatten_list = [item for sublist in c_comp for item in sublist]
    return flatten_list


def gererate_all_c_components(num_var, c_size):
    # generate all possible c components
    c_comp = {}
    for iCsize in range(1, c_size+1):
        sets = list(itertools.combinations( sorted(set(range(num_var))), iCsize))
        c_comp[iCsize] = sets
    return c_comp

