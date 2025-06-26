import sys
def read_fasta(seq):
    try:
        with open('sequences.fasta','r') as f:
            lines = f.readlines()
            seq_dict = {}
            header = " " 
            for line in lines:
                line = line.strip()
                validation = all(base in 'ATGC' for base in line)
                if line.startswith('>') and validation==False:
                    header = line
                    seq_dict[header] = " "
                elif validation == True: 
                        seq_dict[header] += line
                else:
                    print("Invalid sequence: header misssing or invalid dna sequence")
                    return
            return seq_dict
    except FileNotFoundError:
        print("File not found")
def filter_sequences(seq):
    filtered = {}
    min_len = int(input("Enter the minimum length of your fasta sequence: "))
    for seq_id, dna in seq.items():
        if len(dna) >= min_len:
            filtered[seq_id] = dna
        return filtered
def write_fasta(seq):
    try:
        with open("filtered.fasta","w") as file:
            for seq_id, dna in seq.items():
                file.write(f"{seq_id}\n")
                file.write(f"{dna}\n")
            print("File created in your folder!")     
    except PermissionError:
        print("Error: You don't have permission to write here.")
def main(fasta_file):
    fasta = read_fasta(fasta_file)
    print("Total sequences read: ",len(fasta))
    filtered = filter_sequences(fasta)
    print("Filtered sequences: ",len(filtered))
    write_fasta(filtered)
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("for using this script: sequences.fasta <fasta sequence>")
        sys.exit(1)
    sequences = sys.argv[1] 
    main(sequences)
                    