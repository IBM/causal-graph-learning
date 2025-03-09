# © Copyright IBM Corporation 2025. All Rights Reserved.
# LICENSE: Eclipse Public License - v 2.0, https://opensource.org/licenses/EPL-2.0
# SPDX-License-Identifier: EPL-2.0


# To run, type:
# python3 generate_ccomponents.py sample_file_name c_size single_parent_size other_c_parent_size data_directory
# For example:
# python3 generate_ccomponents.py sample_graph_random_bowfree_10nodes-0.txt 2 3 1 ../data/test/


from test_latent_scores import generate_scores_bidirect
import numpy as np
import sys




def main():
        data_directory = ''
        datasets = ['']

        c_size = 3
        single_parent_size = 1
        other_c_parent_size = 1

        if len(sys.argv) > 1:
                datasets = [sys.argv[1]]
        if len(sys.argv) > 2:
                c_size = int(sys.argv[2])
        if len(sys.argv) > 3:
                single_parent_size = int(sys.argv[3])
        if len(sys.argv) > 4:
                other_c_parent_size = int(sys.argv[4])
        if len(sys.argv) > 5:
                data_directory = sys.argv[5]


        for dataset in datasets:

	        with open(data_directory + dataset, 'rb') as f:
	                data = np.loadtxt(f, skiprows=0)

	        observed_data = data


	        file_name = data_directory + 'score_' + dataset[:-4] + '-cs' + str(c_size) + '_sps' + str(single_parent_size) + '_ocps' + str(other_c_parent_size) + '.pkl'
	        print(file_name)


	        generate_scores_bidirect(observed_data,
	                                 single_c_parent_size = single_parent_size,
	                                 other_c_parent_size = other_c_parent_size,
	                                 c_size = c_size,
	                                 file_name = file_name)


if __name__ == "__main__":

    main()
