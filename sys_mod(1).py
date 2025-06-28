import sys 
if len(sys.argv) != 3:
    sys.exit()
dna = sys.argv[1]
rna = sys.argv[2]
print('The dna seq is ', dna)
print('The length of dna seq is ', len(dna))
print('The rna seq is ', rna)
print('The length of rna seq is ', len(rna))
    
    