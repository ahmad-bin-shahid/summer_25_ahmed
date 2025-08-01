import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Load your filtered CSV file
df = pd.read_csv("snv.csv")

# Load reference genome into a dictionary
genome = SeqIO.to_dict(SeqIO.parse("GCF_000001405.40_GRCh38.p14_genomic.fna", "fasta"))

# Chromosome name mapping
chrom_map = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
    "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
    "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9",  "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
    "MT": "NC_012920.1"
}

# Extract ±5 bp context (10 bp total)
context_seqs = []
window_size = 5

for index, row in df.iterrows():
    chrom = str(row["Chrom"])
    pos = int(row["Pos"])

    if chrom not in chrom_map:
        context_seqs.append("N" * 10) 
        continue

    ref_chrom = chrom_map[chrom]

    try:
        chrom_seq = genome[ref_chrom].seq
        # Convert from 1-based VCF to 0-based Python indexing
        start = max(0, pos - window_size - 1)
        end = pos + window_size - 1  # 10 bp total (±5 bp)
        context = chrom_seq[start:end]
        context_seqs.append(str(context))
    except KeyError:
        context_seqs.append("N" * 10)

df["ContextSeq"] = context_seqs

# Save to CSV
df.to_csv("mapped_variants.csv", index=False)
print("Saved 10bp window sequences to 'mapped_variants.csv'")
