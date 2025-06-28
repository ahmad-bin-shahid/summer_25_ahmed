import sys
if len(sys.argv) !=2:
    sys.exit()
seq = sys.argv[1]
len_seq = len(seq)
gc_content = (seq.count('G') + seq.count('C')) / len_seq
A_count = seq.count('A')
T_count = seq.count('T')
C_count = seq.count('C')
G_count = seq.count('G')
is_valid = all(base in 'ATGC' for base in seq)
if is_valid:
    print('Sequence is valid')
    print('Sequence length: ', len_seq)
    print('GC content: ', gc_content)
    print('Count of A: ', A_count, '\nCount of T: ', T_count, '\nCount of C: ', C_count, '\nCount of G: ', G_count)
    if gc_content > 0.4:
        print('GC content is higher than threshold (0.4)')
    else:
        print('GC content is lower than threshold (0.4)')
else:
    print('Sequence is not valid')
index = len_seq - 1
reversed_seq = ''
while index>=0:
    reversed_seq += seq[index]
    index -= 1
print("Reversed sequence: ", reversed_seq)