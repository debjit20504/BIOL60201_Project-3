# BIOL60201 Project 3: Python Proteome Analysis Pipeline

## Overview
This repository contains a Python-based proteome analysis pipeline developed as part of the BIOL60201 module. The project involves processing a bacterial genome sequence to perform proteomics-related tasks, including ORF prediction, peptide digestion, mass analysis, and ion statistics calculations. The pipeline automates these tasks through a series of Python scripts executed via a Bash script.

## Repository Structure
- **genome.fasta**: Input nucleotide FASTA file containing the genome sequence of the bacterial organism.
- **group9_task1.py**: Python script for predicting protein ORFs from the genome.
- **group9_task2.py**: Python script for digesting protein sequences into peptides using specified enzymes.
- **group9_task3.py**: Python script for calculating mass-to-charge (m/z) ratios for peptides.
- **group9_task4.py**: Python script for calculating ion statistics and generating tabular outputs.
- **proteome.sh**: Bash script to run the entire pipeline sequentially.

## Pipeline Workflow
The pipeline consists of the following tasks:

### Task 1: ORF Finder
- **Input**: A nucleotide FASTA file (`genome.fasta`).
- **Output**: Predicted protein ORFs in FASTA format.
- **Description**: Identifies all open reading frames (ORFs) in six reading frames, starting with ATG and ending with stop codons (TGA, TAA, TAG). Outputs include frame, length, and start position.

### Task 2: Protein Digester
- **Input**: Protein FASTA file (output from Task 1).
- **Output**: Smaller peptide sequences in a custom format.
- **Description**: Digests protein sequences using user-selected enzymes (Trypsin, Endoproteinase Lys-C, Endoproteinase Arg-C, or V8 proteinase). Supports options for missed cleavages and minimum peptide length.

### Task 3: Mass Analyser
- **Input**: Peptide sequences (output from Task 2).
- **Output**: Table of peptides with mass-to-charge (m/z) values.
- **Description**: Calculates m/z values for peptides, considering user-specified mass types (monoisotopic or average) and charge states.

### Task 4: Ion Statistics Calculator
- **Input**: Table of peptides with m/z values (output from Task 3).
- **Output**: Multiple TSV files with calculated statistics.
- **Description**: Performs statistical analysis on peptide m/z values, including:
  - Counting peptides within specific m/z ranges.
  - Generating binned histograms.
  - Identifying unique peptides.

## How to Run the Pipeline
1. Ensure all scripts and the `genome.fasta` file are in the same directory.
2. Make the Bash script executable:
   ```bash
   chmod +x proteome.sh *.py
   ```
3. Execute the pipeline
   ```bash
   bash proteome.sh
   ```

### Pipeline Options
The Bash script includes randomization for some parameters:
- ORF length threshold.
- Enzyme selection.
- Minimum peptide length and missed cleavages.
- Mass type and charge state.

The parameters are displayed during execution.

### Example Output
The pipeline generates the following output files:
1. task1_output.txt: Predicted protein ORFs.
2. task2_output.txt: Digested peptides.
3. task3_output.txt: Peptide mass-to-charge data.
4. task4_output_*.tsv: Ion statistics in tabular format.

### Requirements
- Python 3.x
- Bash shell
- Standard Python libraries (e.g., argparse, re, panda)

### Credits
This project was completed as part of the BIOL60201 module. The scripts were developed collaboratively by Group 9 members named **Joseph Lee**, **Debjit Pramanik**, **Jawariya** and **Ziyu Yuan**.