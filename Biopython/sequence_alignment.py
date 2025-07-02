#importing relevant libraries 
import sys
from Bio import SeqIO
from Bio import Align
from Bio import pairwise2
from Bio.Seq import Seq
#Function for pairwise alignment:
def align_seq(seq1, seq2, match=1, mis_match=-1, gap_open=-0.5, gap_extended=-0.5):
    align = pairwise2.align.globalms(seq1, seq2, match, mis_match, gap_open, gap_extended)
    best_alignment = align[0] #multiple alignment results in list, getting the best one at first index
    print("Best alignment score is: ", best_alignment.score)
    print("aligned seq 1 is: ", best_alignment.seqA)
    print("aligned seq 2 is: ", best_alignment.seqB)
    print("The start of alignment is: ", best_alignment.start)
    print("The end of alignment is: ", best_alignment.end)
    return best_alignment
#Function for similarity search between 2 sequences:
def similarity(alignment):
    seq_1 = alignment.seqA
    seq_2 = alignment.seqB
    start = alignment.start
    end = alignment.end
    #aligned regions in both sequences
    aligned1 = seq_1[start:end] 
    aligned2 = seq_2[start:end]
    matches = 0
    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i] and aligned1[i] != "-":  #Check if the sequence is aligned and there is no gap (for each base)
            matches += 1 
    length = end - start
    similarity = (matches / length) * 100 if length > 0 else 0 #Formula for similarity
    print("the similarity (%) of alignment is: ", similarity, '%')
    return similarity
#Function for checking gap frequencies
def gap_freq(alignment):
    sq1 = alignment.seqA
    sq2 = alignment.seqB
    gapsA = sq1.count("-")
    gapsB = sq2.count("-")
    frq_sq_1 = (gapsA / len(sq1)) if len(sq1) > 0 else 0
    frq_sq_2 = (gapsB / len(sq2)) if len(sq2) > 0 else 0
    print("Gap freq for seq 1: ", frq_sq_1)
    print("Gap freq for seq 2: ", frq_sq_2)
#Function for checking the conseerved regions in 2 sequences 
def conserved_region(alignment, threshold=20): #any threshold to get minimum length of conserved region
    s1 = alignment.seqA
    s2 = alignment.seqB
    cons_regions = [] #List of conserved regions i.e "match"
    match = ""
    for i in range(len(s1)):
        if s1[i] == s2[i] and s1[i] != '-' and s2[i] != '-':
            match += s1[i] #if the sequence matches at a certain base and there is no gap, sequence added to match up and until a gap appears
        else:
            if len(match) >= threshold: #if length of matched sequence is > threshold, add to the list
                cons_regions.append(match)
            match = ""
    print("Conserved regions (length ≥", threshold, "):")
    if cons_regions:
        print("####".join(cons_regions))
    else:
        print("No conserved region found above threshold.")
    return cons_regions

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py seq1.fna seq2.fna")
        sys.exit(1)

    # Read sequences from FASTA files
    try:
        seq1_file = sys.argv[1]
        seq2_file = sys.argv[2]
        seq1 = str(SeqIO.read(seq1_file, "fasta").seq)
        seq2 = str(SeqIO.read(seq2_file, "fasta").seq)
    except FileNotFoundError:
        print("Error: One or both input files not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA files: {e}")
        sys.exit(1)

    # Perform alignment and analysis
    alignment = align_seq(seq1, seq2)
    similarity(alignment)
    gap_freq(alignment)
    conserved_region(alignment)