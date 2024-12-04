#!/bin/bash

ORF=$(( RANDOM % 26 + 75 ))

echo -e "\nRunning Task 1: This script takes a nucleotide FASTA file as input and outputs predicted protein ORFs in a FASTA format.\n"
echo "ORF selected by the user: 75"
python group9_task1.py genome.fasta 75 task1_output.txt ## TASK 1 COMMAND

# Randomly selecting ENZYME
Enzymes=("Trypsin" "Endoproteinase Lys-C" "Endoproteinase Arg-C" "V8 proteinase")
random_enzyme_index=$(( RANDOM % ${#Enzymes[@]} )) # generating the random index between 0 and len of enzymes list - 1 (which is 3 in our case)

random_enzyme=${Enzymes[$1]} # selecting the random enzyme

# assigning user enzyme as per the selected enzyme
if [ "$random_enzyme" == "Trypsin" ]; then
    user_enzyme="t"
elif [ "$random_enzyme" == "Endoproteinase Lys-C" ]; then
    user_enzyme="lysc"
elif [ "$random_enzyme" == "Endoproteinase Arg-C" ]; then
    user_enzyme="argc"
elif [ "$random_enzyme" == "V8 proteinase" ]; then
    user_enzyme="v8"
fi

# Randomly selecting minimum peptide length and how many missed clevages user want

min_pep=$(( RANDOM % 6 + 5 ))
missed_clevage=$(( RANDOM % 2))

echo -e "\nRunning Task 2: This script takes protein fasta file and returns a file of smaller peptide sequences\n"
echo "Enzyme selected by the user: $random_enzyme"
echo "Missed clevage: $missed_clevage"
echo "Minimum peptide length: 7"
python group9_task2.py -i task1_output.txt -o task2_output.txt -e $user_enzyme -pm 7 -mc 1 ## TASK 2 COMMAND

# Randomly generating mass type and charge
mass_types=("mono" "aver")
mass_type_idx=$(( RANDOM % ${#mass_types[@]} ))

random_mass_type=${mass_types[mass_type_idx]}

charge=$(( RANDOM % 5 + 1 ))

echo -e "\nRunning Task 3: This script takes peptide sequence file and returns a table with info of m/z of every peptide sequence\n"
echo "Mass type: aver"
echo "Charge: $charge"
python group9_task3.py task2_output.txt -m aver -z 1 -o task3_output.txt ## TASK 3 COMMAND

echo -e "\nRunning Task 4: This script takes a table that have info of m/z of every peptide sequence and returns multiple tsv files\n"
for i in {1..4}; do
    python group9_task4.py task3_output.txt -m $i
    echo
done

mv mode1.tsv ${user_enzyme}_mode1.tsv
mv mode2.tsv ${user_enzyme}_mode2.tsv
mv mode3.tsv ${user_enzyme}_mode3.tsv
mv mode4.tsv ${user_enzyme}_mode4.tsv