from Bio import SeqIO
from Bio.Seq import Seq
#Write genbank file into FASTA file:
seq = SeqIO.parse("sequence.gb","genbank")
with open("gb_to_fasta.fasta",'w') as f:
    for i in seq:
        i.annotations["molecule_type"] = "DNA"
        seq_id = i.id
        seq_features = i.features
        seq_sequence = i
        SeqIO.write(seq_sequence, f, "fasta")
print("EXECUTED")