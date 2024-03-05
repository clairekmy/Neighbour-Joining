#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np


input_file = "hiv_sequences.fasta"
sequences = [record for record in SeqIO.parse(input_file, "fasta")]


max_length = max(len(record.seq) for record in sequences)

# Pad each sequence with "-" 
padded_sequences = [SeqRecord(Seq(str(record.seq).ljust(max_length, '-')), id=record.id, description=record.description) for record in sequences]




def check_data_quality(fasta_file):
    seq_records = list(SeqIO.parse(fasta_file, "fasta"))
    all_gap_sequences = 0
    excessive_padding = 0
    sequence_lengths = []

    for record in seq_records:
        sequence = str(record.seq)
        sequence_lengths.append(len(sequence))

       
        if set(sequence) == {"-"}:
            all_gap_sequences += 1

       
        if sequence.count("-") / len(sequence) > 0.5:
            excessive_padding += 1

    
    if len(set(sequence_lengths)) != 1:
        print("Warning: Inconsistent sequence lengths detected.")

    print(f"Total sequences: {len(seq_records)}")
    print(f"Sequences entirely made up of gaps: {all_gap_sequences}")
    print(f"Sequences with more than 50% gaps: {excessive_padding}")


    
    

check_data_quality("hiv_sequences.fasta")





output_fasta_file = "hiv_sequences_padded.fasta"
SeqIO.write(padded_sequences, output_fasta_file, "fasta")

print(f"Saved the padded sequences as a FASTA file: {output_fasta_file}")


sequences_array = np.array([list(str(seq.seq)) for seq in padded_sequences])

print(sequences_array.shape)


output_npy_file = "hiv_sequences_padded.npy"
np.save(output_npy_file, sequences_array)

print(f"Saved the padded sequences as a NPY file: {output_npy_file}")


# In[2]:


import subprocess
import matplotlib.pyplot as plt
import sys
import numpy as np


dataset_name = 'hiv_sequences_padded.npy'

n_leaves_values = [2000]  


running_times = {'SparseNJ': [], 'FNJ': [], 'RNJ': []}

for n_leaves in n_leaves_values:
    
    result = subprocess.run([sys.executable, 'sparseNJ.py', '--dataset', dataset_name, '--num_of_seed', str(10)], capture_output=True, text=True)
    output = result.stdout.strip().split('\n')
    sparse_nj_time = None
    for line in output:
        if "SNJ Elapsed time:" in line:
            time_str = line.split('[')[1].split(']')[0]
            times = [float(x) for x in time_str.split()]
            sparse_nj_time = np.mean(times)
            break
    if sparse_nj_time is None:
        print(f"Error: Could not find SNJ elapsed time in the output for n_leaves={n_leaves}")
    else:
        running_times['SparseNJ'].append(sparse_nj_time)

    
    result = subprocess.run([sys.executable, 'FNJ.py', '--dataset', dataset_name], capture_output=True, text=True)
    output = result.stdout.strip().split('\n')
    fnj_time = None
    for line in output:
        if "FNJ clustering time:" in line:
            fnj_time = float(line.split(':')[1].strip())
            break
    if fnj_time is None:
        print(f"Error: Could not find FNJ clustering time in the output for n_leaves={n_leaves}")
    else:
        running_times['FNJ'].append(fnj_time)

    
    result = subprocess.run([sys.executable, 'RNJ.py', '--dataset', dataset_name], capture_output=True, text=True)
    output = result.stdout.strip().split('\n')
    rnj_time = None
    for line in output:
        if "RNJ clustering time=" in line:
            rnj_time = float(line.split('=')[1].strip())
            break
    if rnj_time is None:
        print(f"Error: Could not find RNJ clustering time in the output for n_leaves={n_leaves}")
    else:
        running_times['RNJ'].append(rnj_time)


        

plt.figure(figsize=(12, 6))  
x = np.arange(len(n_leaves_values))  
width = 0.15  


offset = width + 0.05  

plt.bar(x - offset, running_times['SparseNJ'], width, label='SparseNJ', color='#2ca02c', edgecolor='black')
plt.bar(x, running_times['FNJ'], width, label='FNJ', color='#d62728', edgecolor='black')
plt.bar(x + offset, running_times['RNJ'], width, label='RNJ', color='#1f77b4', edgecolor='black')


plt.xlabel('Algorithm', fontsize=16, labelpad=12)  
plt.ylabel('Elapsed Time (seconds)', fontsize=16, labelpad=12)  
plt.title('Elapsed Time Comparison', fontsize=18, pad=22)  


plt.yticks(fontsize=14)


plt.legend(fontsize=14, frameon=True, shadow=True)


plt.savefig('elapsed_time_comparison_scaled.png', dpi=300, bbox_inches='tight')


plt.show()
       
        
        
        


# In[3]:


import os
import pickle
import networkx as nx


class Args:
    dataset = 'hiv_sequences_padded.npy'
    num_of_seed = 10  

args = Args()




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


# Directory for pickle files and Newick files
pickle_dir = 'results'
newick_dir = 'results'

# Ensure the Newick directory exists
if not os.path.exists(newick_dir):
    os.makedirs(newick_dir)


num_of_seed =10

for seed in range(num_of_seed):
    pickle_file = os.path.join(pickle_dir, os.path.basename(args.dataset) + f'_SNJ_tree_seed{seed}.pickle')
    newick_file = os.path.join(newick_dir, os.path.basename(args.dataset) + f'_SNJ_tree_seed{seed}.newick')

    
    with open(pickle_file, 'rb') as f:
        tree = pickle.load(f)

    
    newick_tree = graph_to_newick(tree)

    
    with open(newick_file, 'w') as f:
        f.write(newick_tree)

print('Conversion to Newick format completed.')


# In[4]:


def tree_to_newick(node):
    """Convert a tree to Newick format."""
    if isinstance(node, str):
        return node
    else:
        return '(' + ','.join(tree_to_newick(child) for child in node) + ')'


with open('results/RNJ_tree.pickle', 'rb') as file:
    rnj_tree = pickle.load(file)


newick_tree = tree_to_newick(rnj_tree) + ';'


with open('results/RNJ_tree.newick', 'w') as file:
    file.write(newick_tree)


# In[5]:


def tree_to_newick(node):
    """Convert a tree to Newick format."""
    if isinstance(node, str):
        return node
    else:
        return '(' + ','.join(tree_to_newick(child) for child in node) + ')'


with open('results/FNJ_tree.pickle', 'rb') as file:
    fnj_tree = pickle.load(file)


newick_tree = tree_to_newick(fnj_tree) + ';'


with open('results/FNJ_tree.newick', 'w') as file:
    file.write(newick_tree)


# In[6]:


import matplotlib.pyplot as plt
import dendropy
import numpy as np


treefile_path = "hiv_sequences_padded.fasta.treefile"
reference_tree = dendropy.Tree.get(
    path=treefile_path,
    schema="newick"
)



rnj_tree_path = "results/RNJ_tree.newick"
rnj_tree = dendropy.Tree.get(
    path=rnj_tree_path,
    schema="newick"
)

fnj_tree_path = "results/FNJ_tree.newick" 
fnj_tree = dendropy.Tree.get(
    path=fnj_tree_path,
    schema="newick"
)

snj_trees = []
num_of_seeds = 10  
for seed in range(num_of_seeds):
    snj_tree_path = f"results/hiv_sequences_padded.npy_SNJ_tree_seed{seed}.newick"
    snj_tree = dendropy.Tree.get(
        path=snj_tree_path,
        schema="newick"
    )
    snj_trees.append(snj_tree)


taxon_namespace = dendropy.TaxonNamespace()


reference_tree.migrate_taxon_namespace(taxon_namespace)
rnj_tree.migrate_taxon_namespace(taxon_namespace)
fnj_tree.migrate_taxon_namespace(taxon_namespace)
for snj_tree in snj_trees:
    snj_tree.migrate_taxon_namespace(taxon_namespace)


rnj_transfer_distance = dendropy.calculate.treecompare.symmetric_difference(reference_tree, rnj_tree)
fnj_transfer_distance = dendropy.calculate.treecompare.symmetric_difference(reference_tree, fnj_tree)
snj_transfer_distances = [
    dendropy.calculate.treecompare.symmetric_difference(reference_tree, snj_tree)
    for snj_tree in snj_trees
]


mean_snj_transfer_distance = np.mean(snj_transfer_distances)


print(f"Transfer distance (RNJ vs. Reference): {rnj_transfer_distance}")
print(f"Transfer distance (FNJ vs. Reference): {fnj_transfer_distance}")
print(f"Mean transfer distance (SNJ vs. Reference): {mean_snj_transfer_distance}")


# In[7]:





distances = [ mean_snj_transfer_distance, rnj_transfer_distance, fnj_transfer_distance]
labels = ['SNJ', 'RNJ', 'FNJ']


plt.bar(labels, distances, color=['blue', 'green', 'red'])


plt.ylabel('Transfer Distance')
plt.title('Transfer Distance Comparison')


plt.savefig('transfer_distance_comparison.png')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




