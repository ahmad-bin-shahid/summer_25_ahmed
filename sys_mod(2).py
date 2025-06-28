import sys
'''
def seq_concat(a,b):
    seq1 = a
    seq2 = b
    concat = seq1 + seq2
    return concat
'''
def gc_content(sequence):
    seq_len = len(sequence)
    GC_Content = ((sequence.count('G') + sequence.count('C')) / seq_len) * 100
    print("The gc content is: ", GC_Content)
if __name__ == "__main__":
    if len(sys.argv) !=2:
        sys.exit("Usage: python file_name arg_1, arg_2")
    seq_1 = sys.argv[1]
    gc_content(seq_1)