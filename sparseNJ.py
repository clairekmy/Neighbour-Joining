import os
import sys
import pickle

sys.path.append(os.getcwd())
import time
import numpy as np
import networkx as nx
from collections import defaultdict
from skbio import DistanceMatrix
from skbio.tree import nj
import copy
import argparse
from utils import convert_data_str_to_onehot, newick2nx


def DFS(src, visited, subtree_size, subtree_leaves, parent_array, n, tree):
    visited[src] = True
    n[0] += 1
    subtree_size[src] = 1

    if src < n_leaves:
        subtree_leaves[src].append(src)

    for adj in tree.adj[src]:
        if not visited[adj] and not centroidMarked[adj]:
            DFS(adj, visited, subtree_size, subtree_leaves, parent_array, n, tree)
            subtree_size[src] += subtree_size[adj]
            parent_array[adj] = src
            subtree_leaves[src] += subtree_leaves[adj]


def getCentroid(src, visited, subtree_size, n, tree):
    is_centroid = True
    visited[src] = True
    heaviest_child = 0

    for adj in tree.adj[src]:
        if not visited[adj] and not centroidMarked[adj]:
            if subtree_size[adj] > n / 2:
                is_centroid = False

            if heaviest_child == 0 or subtree_size[adj] > subtree_size[heaviest_child]:
                heaviest_child = adj

    if is_centroid and n - subtree_size[src] <= n / 2:
        return src

    return getCentroid(heaviest_child, visited, subtree_size, n, tree)


def getCentroidTree(src, tree, subtree_size=None):
    if subtree_size == None:
        visited = [False] * MAXN
        subtree_size = [0] * MAXN
        parent_array = [-1] * MAXN
        subtree_leaves = defaultdict(list)
        n = [0]

        DFS(src, visited, subtree_size, subtree_leaves, parent_array, n, tree)
    else:
        n = [subtree_size[src]]

    visited = [False] * MAXN
    centroid = getCentroid(src, visited, subtree_size, n[0], tree)
    centroidMarked[centroid] = True

    return centroid


def orient_pick(held_out, leaves_at_adj, num_of_ort):
    num_leaf_at_adj = len(leaves_at_adj)
    if num_leaf_at_adj < num_of_ort:
        num_leaf_to_sample = num_of_ort - num_leaf_at_adj
        ort_leaves = list(leaves_at_adj) + list(np.random.choice(leaves_at_adj, num_leaf_to_sample))

    else:
        upp_limit = int(np.log2(n_leaves))
        if num_leaf_at_adj > upp_limit:
            leaves_at_adj_samp = np.random.choice(leaves_at_adj, upp_limit)
        else:
            leaves_at_adj_samp = leaves_at_adj.copy()

        data_leaves_onehot = data_onehot[leaves_at_adj_samp, :, :]
        data_heldout_onehot = np.squeeze(data_onehot[held_out, :, :])

        dist_vector = np.abs(data_leaves_onehot[:, :, :] - data_heldout_onehot[:, :]).sum(axis=-1).sum(axis=-1)
        dist_vector = dist_vector / (2 * data.shape[1])
        dist_vector = (-3 / 4) * np.log(1 - dist_vector * 4 / 3)

        ort_leaves = leaves_at_adj_samp[np.argsort(dist_vector)[0: num_of_ort]]
        ort_leaves = list(ort_leaves)
    return ort_leaves


def decomposeTree(root, tree, held_out, is_last=False, current_degree=0, first_placement=False, subtree_leaves=None,
                  subtree_size=None, parent_array=None, root_parent=-1, dec_thres=0):
    # print('root', root)
    if (
            root, current_degree, tuple(stop_node_dict[root]),
            root_parent) not in root_centroid_dict and current_degree == 0:
        cend_tree = getCentroidTree(root, tree)
    elif (
            root, current_degree, tuple(stop_node_dict[root]),
            root_parent) not in root_centroid_dict and current_degree != 0:
        cend_tree = getCentroidTree(root, tree, subtree_size)
        root_centroid_dict[(root, current_degree, tuple(stop_node_dict[root]), root_parent)] = cend_tree
    elif (root, current_degree, tuple(stop_node_dict[root]), root_parent) in root_centroid_dict:
        cend_tree = root_centroid_dict[(root, current_degree, tuple(stop_node_dict[root]), root_parent)]
        # check if the centroid is still up
        if current_degree != 0:
            curr_adj_sizes = []
            curr_adj_list = []
            for adj in tree.adj[cend_tree]:
                curr_adj_sizes.append(subtree_size[adj])
                curr_adj_list.append(adj)

            sorted_adj_sizes = np.sort(np.array(curr_adj_sizes.copy()))
            # print(sorted_adj_sizes)
            sorted_adj_list = np.array(curr_adj_list.copy())
            sorted_adj_list = sorted_adj_list[np.argsort(np.array(curr_adj_sizes))]
            sorted_adj_sizes[-1] = subtree_size[root] - sorted_adj_sizes[0] - sorted_adj_sizes[1] - 1

            if (sorted_adj_sizes > subtree_size[root] / 2).any():
                cend_tree = sorted_adj_list[np.argmax(sorted_adj_sizes)]
                root_centroid_dict[(root, current_degree, tuple(stop_node_dict[root]), root_parent)] = cend_tree

        centroidMarked[cend_tree] = True

    if first_placement and current_degree == 0:
        visited = [False] * MAXN
        subtree_size = [0] * MAXN
        n = [0]
        parent_array = [-1] * MAXN
        subtree_leaves = defaultdict(list)

        DFS(cend_tree, visited, subtree_size, subtree_leaves, parent_array, n, tree)
        root_centroid_dict[(root, current_degree, tuple(stop_node_dict[root]), root_parent)] = cend_tree
        global first_dfs_subtree_leaves
        global first_dfs_subtree_sizes
        global first_parent_array
        first_dfs_subtree_leaves = copy.deepcopy(subtree_leaves)
        first_dfs_subtree_sizes = copy.deepcopy(subtree_size)
        first_parent_array = copy.deepcopy(parent_array)

    orient_leaves_1 = []
    orient_leaves_2 = []
    orient_leaves_3 = []

    subtree_leaf_size = []
    active_adj = []
    if current_degree == 0:
        subtree_all_leaves = subtree_leaves[cend_tree]
    else:
        subtree_all_leaves = subtree_leaves[root]

    for adj in tree.adj[cend_tree]:
        # print('adj ' + str(adj))
        if adj != parent_array[cend_tree]:
            leaves_at_adj = np.array(subtree_leaves[adj])

            subtree_leaf_size.append(leaves_at_adj.shape[0])
            subtree_all_leaves = list([leaf for leaf in subtree_all_leaves if leaf not in leaves_at_adj])

            if leaves_at_adj.shape[0] == 0:
                leaves_at_adj = np.array(first_dfs_subtree_leaves[adj].copy())

            ort_leaves = orient_pick(held_out, leaves_at_adj, 3)
            orient_leaves_1.append(ort_leaves[0])
            orient_leaves_2.append(ort_leaves[1])
            orient_leaves_3.append(ort_leaves[2])
            active_adj.append(adj)

        else:
            parent_adj = adj

    if len(subtree_all_leaves) > 0:
        subtree_leaf_size.append(len(subtree_all_leaves))
        ort_leaves = orient_pick(held_out, np.array(subtree_all_leaves), 3)
        orient_leaves_1.append(ort_leaves[0])
        orient_leaves_2.append(ort_leaves[1])
        orient_leaves_3.append(ort_leaves[2])
        active_adj.append(parent_adj)
    elif current_degree != 0:
        subtree_leaf_size.append(len(subtree_all_leaves))
        grand_parent = parent_array[root]
        leaves_at_grand_parent = first_dfs_subtree_leaves[grand_parent].copy()
        leaves_at_sibling = list(
            [leaf for leaf in leaves_at_grand_parent if leaf not in first_dfs_subtree_leaves[root]])
        ort_leaves = orient_pick(held_out, np.array(leaves_at_sibling), 3)
        orient_leaves_1.append(ort_leaves[0])
        orient_leaves_2.append(ort_leaves[1])
        orient_leaves_3.append(ort_leaves[2])
        active_adj.append(parent_adj)

    orient_leaves_1.append(held_out)
    orient_leaves_2.append(held_out)
    orient_leaves_3.append(held_out)

    dist_matrix = np.zeros([3, len(orient_leaves_1), len(orient_leaves_1)])

    data_quart1_onehot = data_onehot[orient_leaves_1, :]
    data_quart2_onehot = data_onehot[orient_leaves_2, :]
    data_quart3_onehot = data_onehot[orient_leaves_3, :]
    for ki in range(len(orient_leaves_1)-1):
        dist_matrix[0, ki, ki+1:] = np.abs(data_quart1_onehot[ki+1:, :, :] - data_quart1_onehot[ki, :, :]).sum(axis=-1).sum(
            axis=-1)
        dist_matrix[1, ki, ki+1:] = np.abs(data_quart2_onehot[ki+1:, :, :] - data_quart2_onehot[ki, :, :]).sum(axis=-1).sum(
            axis=-1)
        dist_matrix[2, ki, ki+1:] = np.abs(data_quart3_onehot[ki+1:, :, :] - data_quart3_onehot[ki, :, :]).sum(axis=-1).sum(
            axis=-1)


    dist_matrix = dist_matrix / (2 * data.shape[1])
    dist_matrix = (-3 / 4) * np.log(1 - (dist_matrix * 4 / 3))
    dist_matrix = dist_matrix.mean(axis=0)
    dist_matrix += dist_matrix.T

    dm = DistanceMatrix(dist_matrix)
    NJ_tree = nj(dm)

    held_out_idx = len(orient_leaves_1) - 1
    selected_path = int(NJ_tree.find(str(held_out_idx)).siblings()[0].name)
    selected_adj = active_adj[selected_path]

    if not centroidMarked[selected_adj]:
        if selected_adj == parent_array[cend_tree]:
            grand_parent = parent_array[cend_tree]
            while grand_parent != -1 and grand_parent != parent_array[root]:
                subtree_leaves[grand_parent] = list(
                    [leaf for leaf in subtree_leaves[grand_parent] if leaf not in subtree_leaves[cend_tree]])
                subtree_size[grand_parent] -= subtree_size[cend_tree]
                stop_node_dict[grand_parent].append(cend_tree)
                grand_parent = parent_array[grand_parent]

            subtree_leaves[cend_tree] = []
            subtree_size[cend_tree] = 0

            selected_adj = root

    if selected_adj != root:
        root_parent = cend_tree
    if current_degree == 0:
        root_parent = cend_tree

    current_degree += 1

    if selected_adj >= n_leaves and not centroidMarked[selected_adj]:
        return decomposeTree(selected_adj, tree, held_out, current_degree=current_degree, subtree_leaves=subtree_leaves,
                             subtree_size=subtree_size, parent_array=parent_array, root_parent=root_parent,
                             dec_thres=dec_thres)
    else:
        if selected_adj < n_leaves:
            next_selected_adj = parent_array[selected_adj]
        else:
            next_selected_adj = cend_tree

        return selected_adj, next_selected_adj, current_degree




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run sparse Neighbor Joining on a dataset.')
    parser.add_argument('--dataset', required=True, help='Path to the dataset file (e.g., data/HIV.npy)')
    parser.add_argument('--num_of_seed', type=int, default=10, help='Number of seeds (e.g., 5, 10, 20)')
    args = parser.parse_args()

    # Load the dataset
    data = np.load(args.dataset)

    # Perform analysis here
    print(f"Dataset {args.dataset} loaded with shape {data.shape}")



    data = convert_data_str_to_onehot(data)
    n_leaves = data.shape[0]
    initial_n_leaves = int(np.floor(np.sqrt(n_leaves * np.log2(n_leaves))))
    num_of_heldout = n_leaves - initial_n_leaves
    num_of_seed = args.num_of_seed
    dist_in_edge = np.zeros([num_of_seed, num_of_heldout])
    time_array = np.zeros([num_of_seed])
    seq_length = data.shape[1]
    data_onehot = data

    # Ensure results directory exists
    results_dir = 'results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    for seed in range(num_of_seed):
        np.random.seed(seed)



        MAXN = 2 * n_leaves - 2

        first_dfs_subtree_sizes = [0] * MAXN
        first_dfs_subtree_leaves = defaultdict(list)
        first_parent_array = [-1] * MAXN

        root_centroid_dict = {}
        parent_node_dict = defaultdict()
        parent_adj_dict = defaultdict()

        perm_idx = np.random.permutation(n_leaves)
        first_leaf_at_rand = perm_idx[0]

        # RANDOM OPTION
        initial_leaves = perm_idx[0:initial_n_leaves]
        held_out_leaves = perm_idx[initial_n_leaves:]

        our_tic = time.process_time()
        data_onehot_initial = data_onehot[initial_leaves, :, :]

        dist_matrix = np.zeros([initial_n_leaves, initial_n_leaves])
        for ki in range(initial_n_leaves):
            dist_matrix[ki, ki+1:] = np.abs(data_onehot_initial[ki+1:, :, :] - data_onehot_initial[ki, :, :]).sum(axis=-1).sum(axis=-1)

        dist_matrix += dist_matrix.T
        dist_matrix = dist_matrix / (2 * data.shape[1])
        dist_matrix = (-3 / 4) * np.log(np.clip(1 - dist_matrix * 4 / 3, 1e-10, None))  # Avoid division by zero and negative values

        # Ensure the distance matrix is symmetric and does not contain NaNs
        dist_matrix = (dist_matrix + dist_matrix.T) / 2
        dist_matrix = np.nan_to_num(dist_matrix)

        dm = DistanceMatrix(dist_matrix)
        tree_newick = nj(dm, result_constructor=str)
        tree = newick2nx(tree_newick, initial_n_leaves)



        our_toc = time.process_time()
        output_path = os.path.join(results_dir, os.path.basename(args.dataset) + f'_SNJ_tree_seed{seed}.pickle')
        with open(output_path, 'wb') as f:
            pickle.dump(tree, f)

        time_array[seed] = our_toc - our_tic

    print('*******************')
    print('DONE=>>>>> ' + args.dataset)
    print('SparseNJ computation time:', time_array)
    print('*******************')
