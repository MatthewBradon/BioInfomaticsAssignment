import sys
import os

dna_seq_file = "DnaSeq.fasta"

compliment_nucleotide = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

CodonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def read_dna_seq(file):
    with open(file, "r") as f:
        # Skip the first line since it is the header
        f.readline()
        dna_seq = f.read().replace("\n", "")  



        return dna_seq
    
    return

def compliment_strand(dna_seq):
    return "".join(compliment_nucleotide[n] for n in dna_seq[::-1])


def read_codon(dna_seq, reading_frame_no=0):
    if reading_frame_no not in [0, 1, 2]:
        print("Invalid reading frame")
        return
    
    codon = []
    for i in range(reading_frame_no, len(dna_seq), 3):
        codon.append(dna_seq[i:i+3])
    return codon

def condon_to_amino(codon):
    if codon not in CodonTable:
        return None
    
    return CodonTable[codon]

def get_possible_orfs(aminos):
    start_amino = "M"
    stop_amino = "*"
    orfs = []
    
    # Get all the start and stop positions
    start_positions = [i for i, a in enumerate(aminos) if a == start_amino]
    stop_positions = [i for i, a in enumerate(aminos) if a == stop_amino]
    
    for start in start_positions:
        for stop in stop_positions:
            if stop > start:
                chain = aminos[start:stop + 1]
                orfs.append((start, stop, "".join(aminos[start:stop+1])))
                break  # Stop at the first stop codon
    
    return orfs










def main(dna_seq_file):
    # dna_seq is 5' to 3'
    dna_seq = read_dna_seq(dna_seq_file)
    # compliment_dna_seq is 3' to 5'
    compliment_dna_seq = compliment_strand(dna_seq)


    log_file = open("log.txt", "w")

    log_file.write(f"\nPrimary DNA Sequence (5' to 3'):\n")
    log_file.write(dna_seq + "\n\n")

    for frame in [0,1,2]:
        codons = read_codon(dna_seq, frame)
        amino_acids = [CodonTable[c] for c in codons if c in CodonTable]
        orfs = get_possible_orfs(amino_acids)

        log_file.write(f"Reading Frame {frame+1}\n")
        log_file.write("Codons:\n" + " ".join(codons) + "\n\n")
        log_file.write("Amino Acid Sequence:\n" + "".join(amino_acids) + "\n\n")
        
        log_file.write("Open Reading Frames:\n")
        for start, stop, seq in orfs:
            log_file.write(f"Start: {start*3+1}, Stop: {stop*3+3}, Length: {stop-start+1} aa\n")
            log_file.write(f"AA Sequence: {seq}\n\n")



    log_file.write(f"\n\n\n\n\n\nComplimentary DNA Sequence (3' to 5'):\n")
    log_file.write(compliment_dna_seq + "\n\n")

    # Compilment
    for frame in [0,1,2]:
        codons = read_codon(compliment_dna_seq, frame)
        amino_acids = [CodonTable[c] for c in codons if c in CodonTable]
        orfs = get_possible_orfs(amino_acids)


        log_file.write(f"Reading Frame -{frame+1}\n")
        log_file.write("Codons:\n" + " ".join(codons) + "\n\n")
        log_file.write("Amino Acid Sequence:\n" + "".join(amino_acids) + "\n\n")
        
        log_file.write("Open Reading Frames:\n")
        for start, stop, seq in orfs:
            log_file.write(f"Start: {start*3+1}, Stop: {stop*3+3}, Length: {stop-start+1} aa\n")
            log_file.write(f"AA Sequence: {seq}\n\n")
    
    log_file.close()


if __name__ == "__main__":
    # Take 1 argument as input
    if len(sys.argv) != 2:
        print("Invalid arguments please run the program: python assignment.py <dna_seq_file>")
        sys.exit(1)

    dna_seq_file = str(sys.argv[1])

    # Check if the file exists
    if not os.path.exists(dna_seq_file):
        print("File does not exist")
        sys.exit(1)

    # Check if the file is not a .fasta file
    if not (dna_seq_file.endswith(".fasta") or dna_seq_file.endswith(".fa")):
        print("File is not a .fasta file")
        sys.exit(1)

    

    main(dna_seq_file)



    # if orf is less than 10 in length, it is not considered as a valid orf
    # codone bias
    # BLAST the ORF
