import argparse
import re

# Argparse arguments: Reading in the fasta file. Specify length of the ORF. Output file specification 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_fasta", type=argparse.FileType("r"), help="Reads in FASTA file. Example input = python Task1.py file.fasta 100 file_output.txt")
    parser.add_argument("ORF_length", type=int, help="Please Specify the minimum number of AMINO ACIDS in your ORF when using command line. Using this program all codons containing an unknown codon (N) will be ignored")
    parser.add_argument("output", type=str, help="Output File Name")

    args = parser.parse_args()

    sequences, names = fastaread(args.in_fasta)

    output = open(args.output, "w") if args.output else sys.stdout

    for name, seq in sequences.items():
        orfs = locate_orfs(seq, args.ORF_length)
        found = False
        for orf_seq, start_index, orf_length, orf_frame in orfs:
            protein = translate(orf_seq)
            if "X" not in protein:  # Exclude ORFs containing 'X'
                found = True
                orf_name = f">{name}_{orf_frame}_{start_index + 1:04d}"
                output.write(f"{orf_name}\n")
                output.write(f"{protein}\n")
        if not found:
            print(f"No ORF(s) found for sequence {name}")

    if args.output:
        output.close()


def fastaread(FASTA_file):
    sequences = {}
    names = []
    name = None
    for line in FASTA_file:
        line = line.rstrip("\n")
        if line.startswith(">"):
            header = line[1:]
            name = header.split()[0]
            names.append(name)
            sequences[name] = ""
            print("FASTA file read in ")
        elif name:
            sequences[name] += line.upper()  # Made upper to match FASTA format to the translation table
        else:
            raise ValueError("FASTA file in wrong format")
    return sequences, names


def reverse_complement(seq):
    complement = str.maketrans("ATGCatgcn", "TACGtacgn") #Used instead of dictionary 
    return seq.translate(complement)[::-1]


def reading_frames(seq):
    F1 = seq
    F2 = seq[1:]
    F3 = seq[2:]
    reverse = reverse_complement(seq)
    R1 = reverse
    R2 = reverse[1:]
    R3 = reverse[2:]
    return [F1, F2, F3, R1, R2, R3]


def translate(seq):
    encode = {
        "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA": "T", "ACC": "T",
        "ACG": "T", "ACT": "T", "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
        "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I", "CAA": "Q", "CAC": "H",
        "CAG": "Q", "CAT": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", "CTA": "L", "CTC": "L",
        "CTG": "L", "CTT": "L", "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
        "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A", "GGA": "G", "GGC": "G",
        "GGG": "G", "GGT": "G", "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
        "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y", "TCA": "S", "TCC": "S",
        "TCG": "S", "TCT": "S", "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
        "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"
    }
    
    protein = ""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        protein += encode.get(codon, "X")
    return protein


def locate_orfs(seq, min_length=50):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []  # Empty list of ORFs. Adding to the list
    frames = reading_frames(seq)

    for frame_index, frame in enumerate(frames):
        frame_offset = frame_index % 3  # Offset for forward frames
        reverse_frame_check = frame_index >= 3  # Check if it's a reverse frame
        used_ranges = []  # Store already used ORFs here. So they can't be used again. Fixes Overlap

        for start_index in range(0, len(frame) - 2, 3):
            if frame[start_index:start_index + 3] == start_codon:
                for end_index in range(start_index, len(frame) - 2, 3):
                    if frame[end_index:end_index + 3] in stop_codons:
                        orf_seq = frame[start_index:end_index + 3]

                        if len(orf_seq) >= min_length * 3: # Is it start or reverse frame
                            
                            if not reverse_frame_check:
                                sequence_start = start_index + frame_offset  # For forward frames
                            else:
                                sequence_start = len(seq) - (end_index + 3)  # For reverse frames

                            # Check if this ORF overlaps with any existing ORFs
                            if not any(start_index < used_end and end_index + 3 > used_start 
                                       for used_start, used_end in used_ranges):
                                orfs.append(
                                    (
                                        orf_seq,
                                        sequence_start,
                                        len(orf_seq) // 3,
                                        f"F{frame_index + 1}" if frame_index < 3 else f"R{frame_index - 2}",
                                    )
                                )
                                
                                used_ranges.append((start_index, end_index + 3)) # Add this ORF's range to used_ranges before looping through again

                        break

    print(f"A total of {len(orfs)} ORFs were found")
    return orfs


if __name__ == "__main__":
    main()
    
