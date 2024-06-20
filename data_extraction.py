from Bio import Entrez, SeqIO
from Bio.Seq import Seq  
import numpy as np
import os
from http.client import IncompleteRead

# Create 'data' directory if it doesn't exist
data_dir = "data"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

# Key in email here
Entrez.email = "xxx@gmail.com"

# Organisms sorted by number of rows
organisms = [
    {
        "name": "E_coli",
        "search_term": '("Escherichia coli"[Organism] AND 500:3000[SLEN])',
        "length_range": (500, 3000)
    },
    {
        "name": "Chicken",
        "search_term": '("Gallus gallus"[Organism] AND 500:4500[SLEN])',
        "length_range": (500, 4500)
    },
    {
        "name": "Candida_albicans",
        "search_term": '("Candida albicans"[Organism] AND 500:5000[SLEN])',
        "length_range": (500, 5000)
    },
    {
        "name": "Mouse",
        "search_term": '("Mus musculus"[Organism] AND 500:1500[SLEN])',
        "length_range": (500, 1500)
    },
    {
        "name": "Tomato",
        "search_term": '("Solanum lycopersicum"[Organism] AND 500:1000[SLEN])',
        "length_range": (500, 1000)
    },
    {
        "name": "Yeast",
        "search_term": '("Saccharomyces cerevisiae"[Organism] AND 200:3000[SLEN])',
        "length_range": (200, 3000)
    },
    {
        "name": "Mycobacterium_tuberculosis",
        "search_term": '("Mycobacterium tuberculosis"[Organism] AND 1000:5000[SLEN])',
        "length_range": (1000, 5000)
    },
    {
        "name": "Bacillus_subtilis",
        "search_term": '("Bacillus subtilis"[Organism] AND 500:4500[SLEN])',
        "length_range": (500, 4500)
    },
    {
        "name": "HIV",
        "search_term": ('("Human immunodeficiency virus 1"[Organism] OR "Human immunodeficiency virus 2"[Organism] OR '
                        '"Human immunodeficiency virus"[Organism] OR "Human immunodeficiency virus 3"[Organism] OR '
                        '"Simian-Human immunodeficiency virus"[Organism] OR HIV[All Fields]) AND '
                        '(1000:5000[SLEN] AND "Homo sapiens"[Organism] AND is_nuccore[filter])'),
        "length_range": (1000, 5000)
    },
    {
        "name": "Zebrafish",
        "search_term": '("Danio rerio"[Organism] AND 500:3000[SLEN])',
        "length_range": (500, 3000)
    },
    {
        "name": "Arabidopsis_thaliana",
        "search_term": '("Arabidopsis thaliana"[Organism] AND 500:5000[SLEN])',
        "length_range": (500, 5000)
    },
    {
        "name": "Salmon",
        "search_term": '("Salmo salar"[Organism] AND 1000:5000[SLEN])',
        "length_range": (1000, 5000)
    },
    {
        "name": "Corn",
        "search_term": '("Zea mays"[Organism] AND 1000:5000[SLEN])',
        "length_range": (1000, 5000)
    },
    {
        "name": "Wheat",
        "search_term": '("Triticum aestivum"[Organism] AND 500:1000[SLEN])',
        "length_range": (500, 1000)
    },
    {
        "name": "Aspergillus",
        "search_term": '("Aspergillus"[Organism] AND 1000:5000[SLEN])',
        "length_range": (1000, 5000)
    },
    {
        "name": "Drosophila_melanogaster",
        "search_term": '("Drosophila melanogaster"[Organism] AND 500:1000[SLEN])',
        "length_range": (500, 1000)
    },
    {
        "name": "Shrimp",
        "search_term": '("Penaeus monodon"[Organism] AND 500:5000[SLEN])',
        "length_range": (500, 5000)
    },
    {
        "name": "Human",
        "search_term": '("Homo sapiens"[Organism] AND 1000:1500[SLEN])',
        "length_range": (1000, 1500)
    },
    {
        "name": "Streptococcus_pneumoniae",
        "search_term": '("Streptococcus pneumoniae"[Organism] AND 500:4500[SLEN])',
        "length_range": (500, 4500)
    },
    {
        "name": "Rice",
        "search_term": '("Oryza sativa"[Organism] AND 500:5000[SLEN])',
        "length_range": (500, 5000)
    }
]

# Categories split into five ranges
categories = {
    "100_to_350": (100, 350),
    "400_to_600": (400, 600),
    "700_to_900": (700, 900),
    "1000_to_1200": (1000, 1200),
    "1300_to_1500": (1300, 1500)
}

def fetch_and_save_sequences(organism, search_term, length_range, rows_range, category):
    """Fetch sequences for a given organism and save them to a FASTA file."""
    search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=rows_range[1])
    search_results = Entrez.read(search_handle)
    search_handle.close()
    total_count = int(search_results['Count'])
    retmax = min(rows_range[1], total_count)  # Ensuring maximum rows

    seq_ids = search_results["IdList"][:retmax]
    if total_count > retmax:
        search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=retmax)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        seq_ids = search_results["IdList"][:retmax]

    # Fetching and saving sequences
    attempts = 0
    max_attempts = 3
    valid_sequences = []
    while attempts < max_attempts:
        try:
            fetch_handle = Entrez.efetch(db="nucleotide", id=seq_ids, rettype="fasta", retmode="text")
            sequences = SeqIO.parse(fetch_handle, "fasta")
            for seq_record in sequences:
                seq_len = len(seq_record.seq)
                if length_range[0] <= seq_len <= length_range[1]:
                    seq_str = str(seq_record.seq).upper()
                    if seq_str.count('-') > seq_len / 2 or seq_str.count('N') > seq_len / 2:
                        continue  # Skip sequences with more than 50% gaps or Ns
                    valid_sequences.append(seq_record)
            fetch_handle.close()
            break
        except IncompleteRead:
            attempts += 1
            if attempts == max_attempts:
                raise
            print(f"Attempt {attempts} failed, retrying...")
    output_file = os.path.join(data_dir, f"{organism}.fasta")
    SeqIO.write(valid_sequences, output_file, "fasta")
    print(f"Saved {len(valid_sequences)} sequences for {organism} in category {category} to {output_file}")
    return valid_sequences

def pad_sequences_and_save_as_npy(valid_sequences, organism, category):
    """Pad sequences, save as both FASTA and NPY."""
    max_len = max(len(seq.seq) for seq in valid_sequences)
    padded_sequences = []
    for seq in valid_sequences:
        padded_seq = seq.seq + '-' * (max_len - len(seq.seq))
        padded_sequences.append(padded_seq)

    padded_fasta_file = os.path.join(data_dir, f"{organism}_padded.fasta")
    with open(padded_fasta_file, 'w') as f:
        for seq, original_seq in zip(padded_sequences, valid_sequences):
            SeqIO.write(SeqIO.SeqRecord(seq=Seq(str(seq)), id=original_seq.id, description=""), f, "fasta")

    npy_file = os.path.join(data_dir, f"{organism}.npy")
    sequences_array = np.array([list(str(seq)) for seq in padded_sequences])
    np.save(npy_file, sequences_array)
    print(f"Padded and saved sequences for {organism} in category {category} in {padded_fasta_file} and {npy_file}")

# Fetch and process sequences for each organism and category
num_categories = len(categories)
for i, (category, rows_range) in enumerate(categories.items()):
    organisms_subset = organisms[i::num_categories]
    for organism_data in organisms_subset:
        key, search_term, length_range = organism_data["name"], organism_data["search_term"], organism_data["length_range"]
        valid_sequences = fetch_and_save_sequences(key, search_term, length_range, rows_range, category)
        if valid_sequences:  # Only pad and save if sequences are found
            pad_sequences_and_save_as_npy(valid_sequences, key, category)

# Check the shape of each dataset
datasets = [
    'data/E_coli.npy', 'data/Chicken.npy', 'data/Candida_albicans.npy', 'data/Mouse.npy', 'data/Tomato.npy',
    'data/Yeast.npy', 'data/Mycobacterium_tuberculosis.npy', 'data/Bacillus_subtilis.npy', 'data/HIV.npy',
    'data/Zebrafish.npy', 'data/Arabidopsis_thaliana.npy', 'data/Salmon.npy', 'data/Corn.npy', 'data/Wheat.npy',
    'data/Aspergillus.npy', 'data/Drosophila_melanogaster.npy', 'data/Shrimp.npy', 'data/Human.npy',
    'data/Streptococcus_pneumoniae.npy', 'data/Rice.npy'
]

for dataset in datasets:
    data = np.load(dataset)
    print(f"{dataset[:-4]}: {data.shape}")
