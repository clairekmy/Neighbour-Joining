import matplotlib.pyplot as plt
import dendropy
import numpy as np
import os
import pickle
from Bio import Phylo

datasets = ['HIV.npy', 'SARS-CoV-2.npy', 'Swine H1N1.npy']
num_of_seeds = 10

# Initialize lists to store transfer distances for all datasets
snj_distances_all = []
rnj_distances_all = []
fnj_distances_all = []

for dataset in datasets:
    # Load the reference tree from an XML file
    xml_tree_path = f"MAFFT_Tree/{dataset.split('.')[0]}.xml"
    # Read the PhyloXML tree
    phyloxml_trees = Phylo.read(xml_tree_path, "phyloxml")
    # Convert to Newick and save
    newick_tree_path = f"MAFFT_Tree/{dataset.split('.')[0]}.fasta.treefile"
    Phylo.write(phyloxml_trees, newick_tree_path, "newick")
    # Load the Newick tree into DendroPy
    reference_tree = dendropy.Tree.get(
        path=newick_tree_path,
        schema="newick"
    )

    # Load the RNJ tree from a pickle file
    rnj_tree_path = f"results/RNJ_tree_data/{dataset}.pickle"
    with open(rnj_tree_path, 'rb') as rnj_file:
        rnj_tree_data = pickle.load(rnj_file)
        print("RNJ Tree Newick String:", rnj_tree_data)
        # Assuming rnj_tree_data is a valid Newick string
        rnj_tree = dendropy.Tree.get(data=rnj_tree_data, schema="newick", taxon_namespace=reference_tree.taxon_namespace)

    # Load the FNJ tree from a pickle file
    fnj_tree_path = f"results/FNJ_tree_data/{dataset}.pickle"  
    with open(fnj_tree_path, 'rb') as fnj_file:
        fnj_tree_data = pickle.load(fnj_file)
        # Assuming fnj_tree_data is a valid Newick string
        fnj_tree = dendropy.Tree.get(data=fnj_tree_data, schema="newick", taxon_namespace=reference_tree.taxon_namespace)

    snj_trees = []
    for seed in range(num_of_seeds):
        snj_tree_path = f"IQTree/{dataset}_SNJ_tree_seed{seed}.newick"  # Modified to load from Newick files
        snj_tree = dendropy.Tree.get(
            path=snj_tree_path,
            schema="newick",
            taxon_namespace=reference_tree.taxon_namespace
        )
        snj_trees.append(snj_tree)

    rnj_transfer_distance = dendropy.calculate.treecompare.symmetric_difference(reference_tree, rnj_tree)
    fnj_transfer_distance = dendropy.calculate.treecompare.symmetric_difference(reference_tree, fnj_tree)
    snj_transfer_distances = [
        dendropy.calculate.treecompare.symmetric_difference(reference_tree, snj_tree)
        for snj_tree in snj_trees
    ]

    mean_snj_transfer_distance = np.mean(snj_transfer_distances)

    # Append the distances for this dataset to the corresponding lists
    snj_distances_all.append(mean_snj_transfer_distance)
    rnj_distances_all.append(rnj_transfer_distance)
    fnj_distances_all.append(fnj_transfer_distance)

# Plotting all datasets in the same graph
labels = [dataset.split('.')[0] for dataset in datasets]
x = np.arange(len(labels))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width, snj_distances_all, width, label='SNJ', color='blue')
rects2 = ax.bar(x, rnj_distances_all, width, label='RNJ', color='green')
rects3 = ax.bar(x + width, fnj_distances_all, width, label='FNJ', color='red')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Transfer Distance')
ax.set_title('Transfer Distance Comparison Across Datasets')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

plt.savefig('transfer_distance_comparison_all_datasets.png')
plt.show()
