# Phylogenetic Tree Construction Using Neighbor Joining Methods

This project focuses on the implementation of various neighbor joining methods for phylogenetic tree construction, including Relaxed Neighbor Joining (RNJ), Fast Neighbor Joining (FNJ), and a novel Sparse Neighbor Joining (SNJ) method. The SNJ method is a modification of the algorithm available at [sparseNJ](https://github.com/kurtsemih/sparseNJ/blob/main/sparseNJ.py).

## Dataset

The HIV datasets used in this project were downloaded from the National Center for Biotechnology Information (NCBI). These datasets contain sequences with lengths ranging from 1000 to 9999 nucleotides.

## Hardware

The models and algorithms in this project are optimized for and run using Apple M1 chips, leveraging their high performance and efficiency for computational tasks.

## Files

- `sparseNJ.py`: Modified version of the Sparse Neighbor Joining algorithm. Ensures the distance matrix is valid and contains no null values.

- `RNJ.py`: Implementation of the Relaxed Neighbor Joining algorithm. Also ensures the distance matrix is valid and contains no null values.

- `FNJ.py`: Implementation of the Fast Neighbor Joining algorithm. Similar to other NJ files, it ensures the distance matrix is valid and contains no null values.

- `final.py`: 
  - Loads the FASTA file and determines the length of the longest sequences, then pads with `-` for the missing gaps.
  - Checks to ensure that there are no sequences made up of entirely gaps and no sequences with more than 50% gaps.
  - Runs the three NJ files and measures the time taken for every model using the same HIV dataset before visualization.
  - Converts all pickle files to Newick format, then measures against a reference tree generated by IQ-TREE before visualization.
 
- 'utils.py'
  - Ensure that the input matrix is converted via one-hot encoding.
  - Converts Newick to NetworkX.
