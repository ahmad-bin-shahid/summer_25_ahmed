import sys
def parse(dna):
    with open('sequences.fasta', 'r') as f:
        try: 
            lines = f.readlines()
            #lines.strip() # removes spaces, additional characters like (\n)
            header = lines[0]
            seq = " "
            for i in lines[1:]: #if you have 1 sequence in fasta file
                seq += i
            print(seq)
        except FileNotFoundError:
            print("File not found")
        return
def read_gff(gff_file):
    with open(gff_file , 'r') as g:
        read = g.readlines()
    with open('neww_genome.gff3', 'w') as wr:
        for i in read:            
            col = i.strip()
            col = i.split('\t')
            id = col[0]
            typ = col[2]
            start_coord = col[3]
            end_coord = col[4]
            score = col[5]
            line = f"{id}\t,{typ}\t,{start_coord}\t,{end_coord}\t,{score}\n"
            wr.write(line)
        print("File created in the same folder")

    
if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('Invalid argument')
    sequence = sys.argv[1]
    read_gff(sequence)
    
        