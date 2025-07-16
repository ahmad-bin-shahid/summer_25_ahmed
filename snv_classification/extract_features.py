import pandas as pd
from Bio.SeqUtils import gc_fraction
from math import log2

# One-hot encode a single base
def one_hot_nuc(nuc):
    return {
        'A': [1, 0, 0, 0],
        'C': [0, 1, 0, 0],
        'G': [0, 0, 1, 0],
        'T': [0, 0, 0, 1]
    }.get(nuc.upper(), [0, 0, 0, 0])  # handle 'N' or invalid bases

# One-hot encode the sequence
def one_hot_sequence_dict(seq):
    seq = seq.upper()
    features = {}
    for i, base in enumerate(seq):
        encoded = one_hot_nuc(base)
        features[f"Base_{i}_A"] = encoded[0]
        features[f"Base_{i}_C"] = encoded[1]
        features[f"Base_{i}_G"] = encoded[2]
        features[f"Base_{i}_T"] = encoded[3]
    return features

# Shannon entropy
def sequence_entropy(seq):
    seq = seq.upper()
    counts = {base: seq.count(base) / len(seq) for base in "ACGT"}
    return -sum(p * log2(p) for p in counts.values() if p > 0)

# Main function
def extract_features(df):
    features_list = []

    for _, row in df.iterrows():
        seq = row['ContextSeq']
        ref = row['Ref']
        alt = row['Alt']

        f = {}

        # One-hot sequence encoding
        f.update(one_hot_sequence_dict(seq))

        # GC/AT content
        f["GC_Content"] = gc_fraction(seq)
        f["AT_Content"] = (seq.count("A") + seq.count("T")) / len(seq)

        # Entropy
        f["Entropy"] = sequence_entropy(seq)

        # Transition/transversion
        transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
        f["Is_Transition"] = int((ref, alt) in transitions)
        f["Label"] = row["Label"]

        features_list.append(f)

    return features_list