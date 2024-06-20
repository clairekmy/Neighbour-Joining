# Neighbor Joining Methods in Phylogenetic Tree Construction

## Overview

This repository contains implementations and comparisons of various Neighbor Joining (NJ) methods used for phylogenetic tree construction. These methods include:

- **Standard Neighbor Joining (NJ)**
- **Relaxed Neighbor Joining (RNJ)**
- **Sparse Neighbor Joining (SparseNJ)**
- **Fast Neighbor Joining (FNJ)**

## Methods

### Neighbor Joining (NJ)
The classic NJ method constructs phylogenetic trees based on distance matrices, minimizing the total branch length.

### Relaxed Neighbor Joining (RNJ)
RNJ offers improved accuracy over NJ by relaxing constraints during tree construction.

### Sparse Neighbor Joining (SparseNJ)
SparseNJ is a modified version of NJ designed to handle large datasets efficiently.

### Fast Neighbor Joining (FNJ)
FNJ accelerates tree construction while maintaining accuracy, making it suitable for large-scale phylogenetic analysis.

## Installation

Clone the repository and install the required packages:
```bash
git clone https://github.com/clairekmy/Neighbour-Joining.git
cd Neighbour-Joining
pip install -r requirements.txt
```

## Usage

Run the different NJ methods using the provided Python scripts. Example:
```bash
python NJ.py
```

## Scripts

- **NJ.py**: Implements the standard Neighbor Joining method.
- **RNJ.py**: Implements the Relaxed Neighbor Joining method.
- **sparseNJ.py**: Implements the Sparse Neighbor Joining method.
- **FNJ.py**: Implements the Fast Neighbor Joining method.
- **data_extraction.py**: Extracts and preprocesses data for analysis.
- **distance.py**: Computes distance matrices used for phylogenetic tree construction.
- **load_distance_matrix.py**: Loads precomputed distance matrices.
- **reference_construction.py**: Constructs reference trees for comparison.
- **time.py**: Measures the execution time of different NJ methods.
- **utils.py**: Utility functions for data manipulation and analysis.

## Data

The repository includes extraction of twenty biological datasets from National Center for Biotechnology Information. These datasets are used to demonstrate the performance of the different NJ methods.

## Results

The performance of each NJ method is compared using metrics such as accuracy and computational time.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to discuss changes.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

For any inquiries, please contact [Claire Mingyue Kang](https://github.com/clairekmy).

