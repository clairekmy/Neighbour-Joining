import matplotlib.pyplot as plt
import dendropy
import numpy as np
import os
import pickle
from skbio import TreeNode
import networkx as nx
from io import StringIO

def graph_to_newick(G, root=None):
    if root is None:
        root_candidates = [n for n, d in G.degree() if d == 2]
        if not root_candidates:
            root_candidates = list(G.nodes())
        if not root_candidates:
            raise ValueError("The graph is empty.")
        root = root_candidates[0]

    def build_newick(node, parent=None):
        children_strs = []
        for neighbor in G.neighbors(node):
            if neighbor == parent:
                continue
            child_str = build_newick(neighbor, node)
            children_strs.append(child_str)
        if children_strs:
            return '(' + ','.join(children_strs) + ')' + str(node)
        else:
            return str(node)

    newick_str = build_newick(root) + ';'
    return newick_str

datasets = [
    'E_coli.npy', 'Yeast.npy', 'Arabidopsis_thaliana.npy', 'Drosophila_melanogaster.npy', 'Salmon.npy',
    'Chicken.npy', 'Mycobacterium_tuberculosis.npy', 'Shrimp.npy', 'Candida_albicans.npy', 'Corn.npy',
    'Bacillus_subtilis.npy', 'Human.npy', 'Mouse.npy', 'HIV.npy', 'Streptococcus_pneumoniae.npy',
    'Wheat.npy', 'Aspergillus.npy', 'Tomato.npy', 'Zebrafish.npy', 'Rice.npy'
]

num_of_seeds = 10

# Initialize lists to store RF distances for all datasets
snj_distances_all = []
rnj_distances_all = []
fnj_distances_all = []
nj_distances_all = []

# Initialize list to store successfully processed datasets
processed_datasets = []

# Initialize lists to store reference trees in Newick format
reference_trees_newick = []

for dataset in datasets:
    # Replace spaces with underscores and remove the file extension for reference tree
    ref_tree_dataset_name = os.path.splitext(dataset)[0].replace(' ', '_')
    
    # Construct the reference tree file path
    ref_tree_path = f"Reference/{ref_tree_dataset_name}.tree"
    if not os.path.exists(ref_tree_path):
        print(f"Reference tree file {ref_tree_path} not found.")
        continue

    reference_tree = dendropy.Tree.get(path=ref_tree_path, schema="newick")

    # Convert reference tree to Newick format and save it
    reference_tree_newick = reference_tree.as_string('newick')
    reference_trees_newick.append(reference_tree_newick)

    # Use original dataset name with spaces for other trees
    dataset_name = os.path.splitext(dataset)[0]
    
    # Load the RNJ tree from a Newick file
    rnj_tree_path = f"results/RNJ_tree_distance_matrix_{dataset_name}.npy.pickle.newick"
    if os.path.exists(rnj_tree_path):
        rnj_tree = dendropy.Tree.get(path=rnj_tree_path, schema="newick", taxon_namespace=reference_tree.taxon_namespace)
    else:
        print(f"RNJ tree file {rnj_tree_path} not found.")
        continue

    # Load the FNJ tree from a Newick file
    fnj_tree_path = f"results/FNJ_tree_distance_matrix_{dataset_name}.npy.pickle.newick"
    if os.path.exists(fnj_tree_path):
        fnj_tree = dendropy.Tree.get(path=fnj_tree_path, schema="newick", taxon_namespace=reference_tree.taxon_namespace)
    else:
        print(f"FNJ tree file {fnj_tree_path} not found.")
        continue

    # Load the NJ tree from a pickle file
    nj_tree_path = f"results/NJ_tree_{dataset_name}.pickle"
    if os.path.exists(nj_tree_path) and os.path.getsize(nj_tree_path) > 0:
        with open(nj_tree_path, 'rb') as pickle_file:
            try:
                nj_tree = pickle.load(pickle_file)
            except (EOFError, pickle.UnpicklingError):
                print(f"Error loading {nj_tree_path}: file is empty or corrupted.")
                continue
    else:
        print(f"NJ tree file {nj_tree_path} not found or empty.")
        continue
    
    # Convert NJ tree to Newick format
    if isinstance(nj_tree, TreeNode):
        newick_buffer = StringIO()
        nj_tree.write(newick_buffer, format='newick')
        nj_newick_str = newick_buffer.getvalue()
    else:
        raise TypeError(f"Expected a scikit-bio TreeNode, but got {type(nj_tree)}")

    nj_tree = dendropy.Tree.get(data=nj_newick_str, schema="newick", taxon_namespace=reference_tree.taxon_namespace)

    # Load SNJ trees from pickle files
    snj_trees = []
    for seed in range(num_of_seeds):
        pickle_tree_path = f"results/{dataset_name}.npy_SNJ_tree_seed{seed}.pickle"
        if os.path.exists(pickle_tree_path) and os.path.getsize(pickle_tree_path) > 0:
            with open(pickle_tree_path, 'rb') as pickle_file:
                try:
                    snj_tree = pickle.load(pickle_file)
                except (EOFError, pickle.UnpicklingError):
                    print(f"Error loading {pickle_tree_path}: file is empty or corrupted.")
                    continue
        
            if isinstance(snj_tree, nx.Graph):
                snj_newick_str = graph_to_newick(snj_tree)
            else:
                raise TypeError(f"Expected a NetworkX graph, but got {type(snj_tree)}")

            snj_tree = dendropy.Tree.get(data=snj_newick_str, schema="newick", taxon_namespace=reference_tree.taxon_namespace)
            snj_trees.append(snj_tree)
        else:
            print(f"SNJ tree file {pickle_tree_path} not found or empty.")
            continue

    # Calculate RF distances
    rnj_rf_distance = dendropy.calculate.treecompare.symmetric_difference(reference_tree, rnj_tree) / 2
    fnj_rf_distance = dendropy.calculate.treecompare.symmetric_difference(reference_tree, fnj_tree) / 2
    nj_rf_distance = dendropy.calculate.treecompare.symmetric_difference(reference_tree, nj_tree) / 2
    snj_rf_distances = [dendropy.calculate.treecompare.symmetric_difference(reference_tree, snj_tree) / 2 for snj_tree in snj_trees]

    mean_snj_rf_distance = np.mean(snj_rf_distances)

    # Get the number of leaves in the reference tree
    num_leaves = len(reference_tree.leaf_nodes())
    max_rf_distance = 2 * (num_leaves - 3)  # for unrooted trees

    # Print debug information
    print(f"Dataset: {dataset}")
    print(f"Number of leaves: {num_leaves}")
    print(f"Max RF distance: {max_rf_distance}")
    print(f"RNJ RF distance: {rnj_rf_distance}")
    print(f"FNJ RF distance: {fnj_rf_distance}")
    print(f"NJ RF distance: {nj_rf_distance}")
    print(f"Mean SNJ RF distance: {mean_snj_rf_distance}")

    # Append the distances for this dataset to the corresponding lists
    snj_distances_all.append(mean_snj_rf_distance)
    rnj_distances_all.append(rnj_rf_distance)
    fnj_distances_all.append(fnj_rf_distance)
    nj_distances_all.append(nj_rf_distance)
    
    # Append the dataset to the processed datasets list
    processed_datasets.append(dataset)

# Convert dataset names for display
processed_dataset_names = [dataset.split('.')[0] for dataset in processed_datasets]

# Create a bar plot
bar_width = 0.2
indices = np.arange(len(processed_datasets))

fig, ax = plt.subplots(figsize=(12,6))

bar1 = ax.bar(indices, nj_distances_all, bar_width, label='NJ', color='purple')
bar2 = ax.bar(indices + bar_width, snj_distances_all, bar_width, label='SparseNJ', color='blue')
bar3 = ax.bar(indices + 2 * bar_width, fnj_distances_all, bar_width, label='FNJ', color='red')
bar4 = ax.bar(indices + 3 * bar_width, rnj_distances_all, bar_width, label='RNJ', color='green')

ax.set_xlabel('Dataset', fontsize=16, labelpad=12)
ax.set_ylabel('RF Distance', fontsize=16, labelpad=12)
ax.set_title('RF Distance Comparison Across Datasets', fontsize=18, pad=22)
ax.set_xticks(range(len(processed_dataset_names)))  # Set x-tick positions
ax.set_xticklabels(processed_dataset_names, rotation=45, ha="right")  # Set and rotate x-tick labels
ax.legend()

plt.tight_layout()
plt.savefig('rf_distance_comparison_bar.png', dpi=300, bbox_inches='tight')
plt.show()
