from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import EcoRI
'''
fasta_seq = SeqIO.read('sequences.fasta','fasta') #for single seequnce in fasta file
seq_id = fasta_seq.id
seq_sequence = fasta_seq.seq
print("ID is: ", seq_id)
print("Sequence is: ", seq_sequence)
multi_seq = SeqIO.parse('sequences.fasta','fasta') #for multiple seq in fasta file (genome)
for i in multi_seq:
    seq_id = i.id
    seq_sequence = i.seq
    print("ID is: ", seq_id)
    #print("Sequence is: ", seq_sequence)print("ID is: ", seq_id)


genbank_seq = SeqIO.parse('sequence.gb','genbank') #for genbank files
for i in genbank_seq:
    seq_id = i.id
    seq_features = i.features
    seq_sequence = i.seq
    subseq = seq_sequence[0:50]
    print("ID is: ", seq_id)
    #print("Feature is: ", seq_features)
    print("Sequence is: ", seq_sequence)
    

#biopython operations:

seq = Seq("ATCGCGCGCTCGATCACaaC")
seq_complement = seq.complement()
print("Compliment: ",seq_complement)
rev_complement = seq.reverse_complement()
print("Reverse complement: ",rev_complement)
rna = seq.transcribe()
print("DNA --> RNA: ",rna)
protein = seq.translate()
print("DNA--->protein(direct): ",protein)
A_count = seq.upper().count("A")
print("Number of 'A': ",A_count)
seq_codon = seq.startswith("TAG")
print("Seq startswith stop codon or not? ", seq_codon)

#Write fasta into genbank file:
seq = SeqIO.read("sequences.fasta","fasta")

seq.annotations["molecule_type"] = "DNA"
with open("New_file",'w') as f:
    SeqIO.write(seq, f, "genbank")
    print("EXECUTED")
  
#Find restriction enzymes
seq = SeqIO.parse("bacteria.fa","fasta")
for i in seq:
    print(i.id)
    sequence = i.seq
    cut_site = EcoRI.search(sequence)
    print("the cut site EcoRI are: ", cut_site)
'''  
#Write fasta into genbank file:
seq = SeqIO.parse("sequences.fasta","fasta")
with open("New_file.gb",'w') as f:
    for i in seq:
        i.annotations["molecule_type"] = "DNA"
        sequeunces = i.seq
        SeqIO.write(i, f, "genbank")
    print("EXECUTED")
