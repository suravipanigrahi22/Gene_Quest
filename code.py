# Install required libraries (Uncomment if running in Colab or first time)
!pip install biopython

from Bio.PDB import PDBParser
from collections import Counter
import os
import requests
import matplotlib.pyplot as plt

# Full names of amino acids
aa_full_names = {
    'ALA': 'Alanine', 'ARG': 'Arginine', 'ASN': 'Asparagine', 'ASP': 'Aspartic Acid',
    'CYS': 'Cysteine', 'GLN': 'Glutamine', 'GLU': 'Glutamic Acid', 'GLY': 'Glycine',
    'HIS': 'Histidine', 'ILE': 'Isoleucine', 'LEU': 'Leucine', 'LYS': 'Lysine',
    'MET': 'Methionine', 'PHE': 'Phenylalanine', 'PRO': 'Proline', 'SER': 'Serine',
    'THR': 'Threonine', 'TRP': 'Tryptophan', 'TYR': 'Tyrosine', 'VAL': 'Valine'
}

# Function to download AlphaFold model
def download_af_model(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        filename = f"AF-{uniprot_id}-F1-model_v4.pdb"
        with open(filename, "wb") as f:
            f.write(response.content)
        print(f" Downloaded structure for {uniprot_id}")
        return filename
    else:
        print(f" No structure found for UniProt ID: {uniprot_id}")
        return None

# Function to get DNA sequence from UniProt
def get_gene_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta?format=fasta&includeIsoform=true"
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = response.text
        lines = fasta_data.strip().split('\n')
        seq = ''.join(lines[1:]).upper()
        return seq
    else:
        print(f" Could not fetch DNA/protein sequence for {uniprot_id}")
        return None




# Ask user for UniProt ID
uniprot_id = input("Enter UniProt ID: ").strip().upper()

# Download and parse structure
pdb_file = download_af_model(uniprot_id)
if pdb_file:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(uniprot_id, pdb_file)

    # Count standard amino acids
    amino_acids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == " ":
                    amino_acids.append(residue.get_resname())

    # Total count
    total_count = len(amino_acids)
    print(f"\n Total amino acids in {uniprot_id}: {total_count}")

    # Top 10 amino acids
    aa_counts = Counter(amino_acids)
    print("\n Top 10 amino acids:")
    for aa, count in aa_counts.most_common(10):
        full_name = aa_full_names.get(aa, aa)
        print(f"{aa} ({full_name}): {count}")

    # Bar plot of top 10 amino acids
    top_aa = aa_counts.most_common(10)
    labels = [aa_full_names.get(aa, aa) for aa, _ in top_aa]
    values = [count for _, count in top_aa]

    plt.figure(figsize=(12, 6))
    bars = plt.bar(labels, values, color='teal')
    plt.title(f"Top 10 Amino Acids in {uniprot_id}")
    plt.xlabel("Amino Acid (Full Name)")
    plt.ylabel("Count")
    plt.xticks(rotation=45)
    plt.grid(True, linestyle='--', alpha=0.5)

    for bar in bars:
        height = bar.get_height()
        plt.annotate(f'{height}', xy=(bar.get_x() + bar.get_width() / 2, height),
                     xytext=(0, 3), textcoords="offset points", ha='center', va='bottom')

    plt.tight_layout()
    plt.show()

# Get DNA/protein sequence and compute nucleotide stats
sequence = get_gene_sequence(uniprot_id)
if sequence:
    print(f"\n Sequence length: {len(sequence)} bases")
    base_counts = Counter(base for base in sequence if base in "ATGC")
    a = base_counts.get('A', 0)
    t = base_counts.get('T', 0)
    g = base_counts.get('G', 0)
    c = base_counts.get('C', 0)
    total = a + t + g + c

    if total > 0:
        at = a + t
        gc = g + c
        print(f"\n Nucleotide content:")
        print(f"A: {a} ({(a/total)*100:.2f}%)")
        print(f"T: {t} ({(t/total)*100:.2f}%)")
        print(f"G: {g} ({(g/total)*100:.2f}%)")
        print(f"C: {c} ({(c/total)*100:.2f}%)")
        print(f"AT content: {at} ({(at/total)*100:.2f}%)")
        print(f"GC content: {gc} ({(gc/total)*100:.2f}%)")
    else:
        print(" The sequence might be a protein (not DNA), as it lacks A/T/G/C.")