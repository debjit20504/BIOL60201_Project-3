
#!/usr/bin/python_

import argparse

# Function to calculate peptide mass:
def pep2mass(seq, mass_type):
    mono = {'A': 71.0371, 'C': 103.0092, 'D': 115.0269, 'E': 129.0426, 'F': 147.0684, 'G': 57.0215, 'H': 137.0589,
            'I': 113.0841, 'K': 128.0950, 'L': 113.0841, 'M': 131.0405, 'N': 114.0429, 'P': 97.0528, 'Q': 128.0586,
            'R': 156.1011, 'S': 87.0320, 'T': 101.0477, 'V': 99.0684, 'W': 186.0793, 'Y': 163.0633, '*': 0.0}
    aver = {'A': 71.08, 'C': 103.14, 'D': 115.09, 'E': 129.12, 'F': 147.18, 'G': 57.05, 'H': 137.14, 'I': 113.16,
            'K': 128.17, 'L': 113.16, 'M': 131.19, 'N': 114.10, 'P': 97.12, 'Q': 128.13, 'R': 156.19, 'S': 87.08,
            'T': 101.10, 'V': 99.13, 'W': 186.21, 'Y': 163.18, '*': 0.0}

    if mass_type not in ["mono", "aver"]:
        raise ValueError(f"Invalid mass type '{mass_type}'. Mass type must be 'mono' or 'aver'.")

    mass_dict = mono if mass_type == "mono" else aver
    water_mass = 18.0106 if mass_type == "mono" else 18.0153

    invalid_amino_acids = [aa for aa in seq if aa not in mass_dict]
    if invalid_amino_acids:
        raise ValueError(f"Invalid amino acids found in sequence: {invalid_amino_acids}")

    peptide_mass = sum(mass_dict[aa] for aa in seq) + water_mass
    return peptide_mass


# Function to read FASTA files:
def fastaread(filename):
    order = []
    seqs = []
    with open(filename, "r") as file:
        prot_name, peptide, missed, sequence, enzyme = None, None, None, [], None
        for line in file:
            line = line.strip()
            if line.startswith(">"):  
                if prot_name and sequence:  
                    full_sequence = "".join(sequence)
                    if full_sequence != "*" * len(full_sequence):  
                         seqs.append((prot_name, peptide, missed, full_sequence, enzyme))
                    else:
                        print(f"Skipping protein '{prot_name}' with only stop codons.")
                header_parts = line[1:].split()
                prot_name = header_parts[0]
                peptide = header_parts[2]
                missed = int(header_parts[3])
                
                if len(header_parts) > 4:  
                    enzyme = " ".join(header_parts[4:])  
                else:
                    raise ValueError(f"Enzyme name missing in header: '{line}'")
                
                order.append(prot_name)
                sequence = [] 
            else:  
                sequence.append(line)

        # Process the last protein
        if prot_name and sequence:
            full_sequence = "".join(sequence)
            if full_sequence != "*" * len(full_sequence):
                seqs.append((prot_name, peptide, missed, full_sequence, enzyme))

    return order, seqs


# Main program
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert peptide sequences to mass-to-charge ratios from a FASTA file. Execution example: python3 group9_task3.py task2_o_v8.txt -m mono -z 1 -o output9.txt")
    parser.add_argument("fileName", help="Input FASTA file")
    parser.add_argument("-m", "--mass_type", choices=["mono", "aver"], default="mono", help="Mass type: mono or aver (default: mono)")
    parser.add_argument("-o", "--output", default="pepmasses.out", help="Output file name (default: pepmasses.out)")
    parser.add_argument("-z", "--charge", type=int, default=1, help="Charge state for mass-to-charge calculation (default: 1)")
    args = parser.parse_args()

    inputfile = args.fileName
    mass_type = args.mass_type
    outputfile = args.output
    charge = args.charge

    if charge <= 0:
        raise ValueError("Charge state must be greater than 0.") 

    # Open the output file
    ofile = open(outputfile, 'w')
    ofile.write("#Prot_name\tpeptide\tmass-to-charge\tz\tp\tsequence\tenzyme\n")

    # Use the fastaread function to read sequences and metadata
    order, seqs = fastaread(inputfile)

    # Write the results for each peptide
    for prot_name, peptide, missed, sequence, enzyme in seqs:
        if sequence: 
            peptide_mass = pep2mass(sequence, mass_type)  
            mass_to_charge = (peptide_mass + (charge * 1.007)) / charge  
            ofile.write(f"{prot_name}\t{peptide}\t{mass_to_charge:.4f}\t{charge}\t{missed}\t{sequence}\t{enzyme}\n")
        else:
            print(f"Skipping protein '{prot_name}' with an empty sequence.")

    # Close the output file
    ofile.close()
    print(f"Results written to {outputfile}")
else:
    print("run as module\n")
    
    
