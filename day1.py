'''
dna = 'ATCG'
dna_int = 2
dna_float = 0.002
print('dtype of dna is', type(dna))
print('dtype of dna_int is', type(dna_int))
print('dtype of dna_float is', type(dna_float))

dna = 'ACTGcGGA'
dna_len = len(dna)
print('length of dna is', dna_len)

dna = input('Enter a DNA sequence: ')
print('Seq of dna is \n',  dna)

dna_1 = 'ATGCGCGTGCA'
dna_2 = 'TCGAGCTCGATCAA'
dna_3 = 'ATGCGCGTGCA'
print('seq1 is equal to seq2?', dna_1 == dna_2)
print('seq1 is not equal to seq2?', dna_1 != dna_2)
print('seq 1 is greater than seq 2?', dna_1 > dna_2)
print('seq 1 is less than seq 2?', dna_1 < dna_2)
print('seq1 and seq2 are equal to seq3?', dna_1 and dna_2 == dna_3)
print('seq1 or seq2 are equal to seq3?', dna_1 or dna_2 == dna_3)

dna = 'ACTCGCGCGCTAGCGAG'
threshold = 10
gc_content = dna.count('G') + dna.count('C')
if gc_content > threshold:
    print('GC content is greater than threshold')
elif gc_content == threshold:
    print('GC content is equal to threshold')
elif gc_content < threshold:
    print('Less than threshold')
else:
    print(gc_content)

for i in range(10):
    print(i)

dna = 'ACTCGCGCGCTAGCGAG'
compliment = ''
for i in dna:
    compliment = i + dna
print(compliment)  
'''