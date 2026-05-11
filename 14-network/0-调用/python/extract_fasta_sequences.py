import sys
import pandas as pd

def read_ids(filename):
    ids = pd.read_csv(filename, sep='\t', header=None, usecols=[0])
    return ids[0].tolist()

def filter_fasta(input_fasta, ids, output_fasta):
    with open(input_fasta, 'r') as fasta_file, open(output_fasta, 'w') as outfile:
        write_sequence = False
        for line in fasta_file:
            if line.startswith('>'):
                sequence_id = line.split('>')[1].strip().split(' ')[0]
                if sequence_id in ids:
                    write_sequence = True
                    outfile.write(line)
                else:
                    write_sequence = False
            elif write_sequence:
                outfile.write(line)

def main(input_ids_file, input_fasta_file, output_fasta_file):
    ids = read_ids(input_ids_file)
    filter_fasta(input_fasta_file, ids, output_fasta_file)
    print("Filtered fasta file has been created.")

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python extract_fasta_sequences.py <meta_file> <fasta_file> <output_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])