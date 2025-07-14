import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from itertools import product

def extract_kmer_frequencies(sequence, k=3):
    """Extract frequencies of all possible k-mers."""
    seq = Seq(sequence)
    kmers = [''.join(p) for p in product('ACGT', repeat=k)]
    kmer_counts = {kmer: 0 for kmer in kmers}
    seq_len = len(seq)
    
    for i in range(seq_len - k + 1):
        kmer = str(seq[i:i+k])
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1

    total_kmers = seq_len - k + 1 if seq_len >= k else 1
    kmer_freqs = {kmer: count / total_kmers for kmer, count in kmer_counts.items()}
    return kmer_freqs

def calculate_gc_at_content(sequence):
    """Calculate GC and AT content using Biopython's gc_fraction."""
    seq = Seq(sequence)
    seq_len = len(seq)
    if seq_len == 0:
        return 0.0, 0.0
    gc_content = gc_fraction(seq)
    at_content = 1.0 - gc_content
    return gc_content, at_content

def tata_box_presence(sequence):
    """Check for TATA box motif variants manually."""
    sequence = str(sequence).upper()
    motifs = ["TATAAAA", "TATAATA", "TATATAA", "TATATTA"]
    for i in range(len(sequence) - 6):
        window = sequence[i:i+7]
        if window in motifs:
            return 1
    return 0

def caat_box_presence(sequence):
    """Check for CAAT box motif."""
    sequence = str(sequence).upper()
    return 1 if "CCAAT" in sequence else 0

def cpg_island_score(sequence):
    """Calculate CpG island score: observed/expected CpG."""
    seq = Seq(sequence)
    seq_len = len(seq)
    if seq_len < 2:
        return 0.0
    cg_count = str(seq).count('CG')
    c_count = str(seq).count('C')
    g_count = str(seq).count('G')
    expected_cpg = (c_count * g_count) / seq_len if seq_len > 0 else 0
    return cg_count / expected_cpg if expected_cpg > 0 else 0

def bre_element_count(sequence):
    """Count BRE-like motifs (manually defined variants)."""
    sequence = str(sequence).upper()
    variants = ["CCACGCC", "CGACGCC", "CCAGCC", "CGAGCC"]
    count = 0
    for i in range(len(sequence) - 6):
        window = sequence[i:i+7]
        if window in variants:
            count += 1
    return count

def one_hot_encode_sequence(sequence, max_len=100):
    """One-hot encode the sequence (A, C, G, T) up to max_len."""
    seq = Seq(sequence)
    encoding = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0],
                'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1]}
    one_hot = []
    sequence = str(seq)[:max_len]
    for base in sequence:
        one_hot.extend(encoding.get(base, [0, 0, 0, 0]))
    if len(sequence) < max_len:
        one_hot.extend([0] * (4 * (max_len - len(sequence))))
    return one_hot

def extract_features(csv_file, output_csv):
    """Extract features from CSV file and save to a new CSV."""
    df = pd.read_csv(csv_file)
    features_list = []
    
    for _, row in df.iterrows():
        sequence = str(row['sequence']).upper()
        gene = row['gene']
        activity = row['activity']
        features = {'gene': gene}
        
        # K-mer frequencies
        kmer_freqs = extract_kmer_frequencies(sequence)
        features.update(kmer_freqs)
        
        # GC & AT content
        gc_content, at_content = calculate_gc_at_content(sequence)
        features['gc_content'] = gc_content
        features['at_content'] = at_content
        
        # Dinucleotide repeats
        features['dinucleotide_repeat_score'] = dinucleotide_repeat_score(sequence)
        
        # Motif presence
        features['tata_box'] = tata_box_presence(sequence)
        features['caat_box'] = caat_box_presence(sequence)
        
        # CpG island
        features['cpg_island_score'] = cpg_island_score(sequence)
        
        # BRE count
        features['bre_element_count'] = bre_element_count(sequence)
        
        # Sequence length
        features['sequence_length'] = len(sequence)
        
        # One-hot encoding
        one_hot = one_hot_encode_sequence(sequence)
        for i, val in enumerate(one_hot):
            features[f'one_hot_{i}'] = val
        
        # Target variable
        features['activity'] = activity
        
        features_list.append(features)
    
    features_df = pd.DataFrame(features_list)
    features_df.to_csv(output_csv, index=False)
    return features_df
