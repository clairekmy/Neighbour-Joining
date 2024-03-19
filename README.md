# Phylogenetic Tree Construction Using Neighbor Joining Methods

This project focuses on the implementation of various neighbor-joining methods for phylogenetic tree construction, including Relaxed Neighbor-Joining (RNJ), Fast Neighbor-Joining (FNJ), and a novel Sparse Neighbor-Joining (SNJ) method. The SNJ method modifies the algorithm available at [sparseNJ](https://github.com/kurtsemih/sparseNJ/blob/main/sparseNJ.py).

## Dataset

The HIV datasets used in this project were downloaded from the National Center for Biotechnology Information (NCBI). These datasets contain sequences with lengths ranging from 1000 to 9999 nucleotides.

## Hardware

The models and algorithms in this project are optimized for and run using Apple M1 chips, leveraging their high performance and efficiency for computational tasks.

## Files

- `sparseNJ.py`: Modified version of the Sparse Neighbor-Joining algorithm. It makes sure that the distance matrix is valid and contains no null values.

- `RNJ.py`: Implementation of the Relaxed Neighbor Joining algorithm. It also ensures the distance matrix is valid and contains no null values.

- `FNJ.py`: Implementation of the Fast Neighbor-Joining algorithm. Similar to other NJ files, it ensures the distance matrix is valid and contains no null values.

- `distance.py`:
  - calculates the distance of the individual trees generated by NJ methods concerning the reference tree.

- `time.py`:
  - measures the time taken for various NJ methods.
  
 
- `utils.py`: 
  - Ensure that the input matrix is converted via one-hot encoding.
  - Converts Newick to NetworkX.
