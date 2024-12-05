import os
import argparse

def parse_arguments():
    enzymes = ['t', 'lysc', 'argc', 'v8']
    parser = argparse.ArgumentParser(description='This script processes a miRNA target dataset performs quality control and outputs a summary with some figures figures.')
    parser.add_argument('--input', '-i', required=True, help='Path to the input file containing protein sequence')
    parser.add_argument('--output', '-o', help='Save the output in a file')
    parser.add_argument('--enzyme', '-e', required=True, choices=enzymes, help='Specify the enzyme from one of these trypsin(t), endoproteinase Lys-C(lysc), endoproteinase Arg-C(argc), v8 proteinase(v8)')
    parser.add_argument('--peptide_length', '-pm', default=7, type=int, help='')
    parser.add_argument('--missed_clevage', '-mc', default=1, type=int, help='')
    
    return parser.parse_args()

def peptide_maker(seq: str, cut_site : str, missed_cleavages: int = 1):
    """
    Splits a protein sequence into peptides based on the enzymes provided with support for missed clevages.
    Args:
        seq (str) - Protein sequence
        cut_site (str) - Clevage site
        missed_clevage (int) - Maximum number of missed clevages allowed
        
    Reuturns:
        list : List of peptides with info about missed clevages
    """
    
    peptides = []
    temp = 0 # starting index of the current peptide
    
    basic_peptides = []  # generating peptides based on the clevage rules
    for i in range(len(seq)):
        if seq[i] in cut_site[0]: # clevage site
            if i+1 >= len(seq) or seq[i+1] != "P": # avoiding cutting if the next amino acid after the clevage site is P
                basic_peptides.append((temp, i+1))
                temp = i+1
    
    if temp < len(seq): # adding the remaining part
        basic_peptides.append((temp, len(seq)))
        
    for idx, (start, end) in enumerate(basic_peptides):
        peptides.append((seq[start:end], 0, cut_site[1]))  # Adding the peptide with 0 missed cleavages

        if missed_cleavages > 0 and idx + 1 < len(basic_peptides): # Adding peptides with 1 missed cleavage
            next_start, next_end = basic_peptides[idx + 1]
            combined_peptide = seq[start:next_end]
            peptides.append((combined_peptide, 1, cut_site[1]))
            
    return peptides
        
def digester_output_maker(header : str, seq_info : set, file : str, peplen : int):
    """
    Writes digested peptide sequences to the output file in a specified format.
    Args:
        header (str): The FASTA header of the protein sequence.
                    It identifies the source protein and will be used to label each peptide.
        seq_info (set): A set of tuples containing:
                        - peptide_seq (str): The peptide sequence after digestion.
                        - missed_clevage (int): The number of missed cleavages in the peptide.
                        - enzyme (str): The enzyme used for digestion.
        file (str): Path to the output file where digested peptides will be written.
        peplen (int): Minimum peptide length. Only peptides with lengths greater than or equal
                    to this value will be included in the output.

    Returns:
        None, The function is just writing peptide sequences and metadata directly to the output file.
    """
    counter = 1
    with open(file, 'a') as f:
        for peptide_seq, missed_clevage, enzyme in seq_info:
            if len(peptide_seq) >= peplen:
                final_header = f"{header} peptide {counter} {missed_clevage} {enzyme}"
                f.write(f"{final_header}\n{peptide_seq}\n")
                counter+=1
    f.close()
            
if __name__ == '__main__':
    
    # dictionary of enzymes with their cleavage sites and names
    enzyme = {
    't' : ['KR', 'Trypsin'],
    'lysc' : ['K', 'Endoproteinase Lys-C'],
    'argc' : ['R', 'Endoproteinase Arg-C'],
    'v8' : ['E', 'V8 proteinase']
    }
    
    args = parse_arguments()
    protein_info = {}
    with open(args.input, 'r') as file: # Read the input protein FASTA file and store headers and sequences
        temp = ''
        for i in file.readlines():
            if ">" in i:
                temp = i.strip()
            else:
                protein_info[temp] = i.strip()
    
    if os.path.exists(args.output): # if the output file exists then remove it
        os.remove(args.output)
        
    for header, seq in protein_info.items():
        peptides = peptide_maker(seq, enzyme[args.enzyme], args.missed_clevage) # Generating peptides for the current protein sequence
        digester_output_maker(header, peptides, args.output, args.peptide_length) # Writing the peptides to the output file, ensuring minimum peptide length is satisfied