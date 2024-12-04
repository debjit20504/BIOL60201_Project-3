#!/usr/bin/python
#
#
#  *********************************** READ ME FIRST ********************************
#  This template script is inciomplete. There are sections missing which you
#  can add. You can fill them in to get it to complete the task. However,
#  if you find thats not an easy thing to do, we suggest you make a new python script,
#  and slowly add bits from this one into it, make sure its working, and then 
#  add some more. That way you can make step by step progress. Have a good look 
#  at this code, and make your own pseudocode plan first, before you start coding
#  your own scripts
#  **********************************************************************************
#

import numpy as np

def readData(file,lower,upper):
    """Read the file and select peptides."""
    dataFile = open(file, "r")
    peps  = {}      # dictionary containing the masses, indexed via a unique identifier

## assumes file format like:   YCL049C    5     774.4059  0 0 GVAMGNVK

    for line in dataFile:
        if not line.startswith('#'):
            (protein, pepnum, mass ) = line.split()[0:3]
            pepseq = line.split()[5]
            #print("data = ",protein, pepnum, mass, pepseq)

            if ( float(mass) < lower or float(mass) > upper ):
                continue

#
# make a unique id for future reference for each peptide (why do you think that might be important?)
#
            pepid = "_".join(list((protein,pepnum,pepseq)))
            peps[pepid] = mass
           
    dataFile.close()
    return (peps) 


import argparse
import csv
import matplotlib.pyplot as plt
import pandas as pd
import os

def save_to_csv(data, output_file, headers):
    """Save the data into CSV file."""
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file,delimiter = "\t")
        writer.writerow(headers)  
        writer.writerows(data)  
    print(f"Results saved to {output_file}")
    

def mode_1_total_count(pepmass, range0, range1, output_file):
    """Count the number of peptides."""
    results = []
    count = len(pepmass)
    results.append([f"{range0:.2f} - {range1:.2f}", count])
    print("Total peptides in range:", count)    
    save_to_csv(results, output_file, headers=["Bin Range (m/z)", "Count"])
    

def mode_2_histogram(pepmass, range0, range1, binsize, numbins, output_file):
    """Generate the binned histogram data."""
    results = []
    #get the  peptide numbers in the range of lower and upper bound, and count the number that peptides fall in the bound
    for i in range(numbins):
        lower = range0 + i * binsize
        upper = lower + binsize
        count = sum(lower <= float(mass) < upper for mass in pepmass.values()) 
        results.append([f"{lower:.2f} - {upper:.2f}", count])
    save_to_csv(results, output_file, headers=["Bin Range (m/z)", "Count"])

    
def mode_3_sliding_window(pepmass, range0, range1, window_size, step_size, output_file):
    """Generate a sliding window of the data."""
    results = []
    current_lower = range0
    #ensure the window size is smaller than the upper range when the window slides
    while current_lower + window_size <= range1: 
        current_upper = current_lower + window_size        
        count = sum(current_lower <= float(mass) < current_upper for mass in pepmass.values()) 
        results.append([f"{current_lower:.2f}-{current_upper:.2f}", count])
        current_lower += step_size
    save_to_csv(results, output_file, headers=["Window Range (m/z)", "Count"])


def mode_4_unique_protein(pepmass, output_file): 
    """Get the identified proteins."""
    unique_proteins = []
    # read the window range when count==1
    with open("mode3.tsv", 'r') as mode3:
        mode3_reader = csv.reader(mode3,delimiter='\t')
        next(mode3_reader)  
        for row in mode3_reader:
            window_range, count = row[0], int(row[1])
            if count == 1:  
                start, end = map(float, window_range.split("-"))
                # find the protein name in the window range and combine them with m/s values
                for prot_name, mass in pepmass.items(): 
                    if start <= float(mass) < end:
                        unique_proteins.append((prot_name, mass))
    #use set() to remove duplicate values
    unique_proteins = list(set(unique_proteins)) 
    print(f"Number of unique proteins: {len(unique_proteins)}")
    save_to_csv(unique_proteins, output_file, headers = ["Prot_name","Mass"])

    
def plot_fixed_histogram(file):
    """Plot mode 2 fixed histogram."""
    data = pd.read_csv("mode2.tsv", sep='\t')
    data["Lower Bound"] = data["Bin Range (m/z)"].str.split(" - ").str[0].astype(float)
    #get the name according to the file name and remove redundant names
    base_name = os.path.splitext(os.path.basename(file))[0]
    if base_name.startswith("Task3_o_"):  
        base_name = base_name.replace("Task3_o_", "")

    plt.bar(data["Lower Bound"], data["Count"], width=binsize, align='edge', color='blue', alpha=0.7)
    plt.title(f"{base_name.capitalize()} Mode 2 Fixed Range Histogram")
    plt.xlabel("Mass-to-Charge Ratio (m/z)")
    plt.ylabel("Number of Peptides")
    plt.savefig("fixed_histogram.png")
    print("Fixed range histogram saved as 'fixed_histogram.png'")
    plt.close()
    
def plot_sliding_window(file):
    """Plot mode 3 Sliding window histogram."""
    data = pd.read_csv("mode3.tsv", sep='\t')
    data["Window Start"] = data["Window Range (m/z)"].str.split("-").str[0].astype(float)
    data["Window End"] = data["Window Range (m/z)"].str.split("-").str[1].astype(float)
    data["Window Center"] = (data["Window Start"] + data["Window End"]) / 2  
    
    base_name = os.path.splitext(os.path.basename(file))[0]
    if base_name.startswith("Task3_o_"):  
        base_name = base_name.replace("Task3_o_", "") 
        
    plt.figure(figsize=(20, 6))
    plt.plot( data["Window Center"], data["Count"], color='blue', linestyle='-', alpha=0.7, label="Sliding Window Counts")
    plt.title(f"{base_name.capitalize()} Sliding Window Histogram")
    plt.xlabel("Mass-to-Charge Ratio (m/z)")
    plt.ylabel("Number of Peptides")
    plt.legend()
    plt.savefig("sliding_window.png")
    print("Mode 3 sliding window plot saved as 'sliding_window.png'")

##
##  Main code starts here
##
parser = argparse.ArgumentParser()
parser.add_argument("fileName", help="name of file containing data")     # a positional argument
parser.add_argument("-s", "--start", help="lower bound of m/z range", default=1000.0, type=float)
parser.add_argument("-e", "--end", help="upper bound of m/z range", default=1500.0, type=float)
parser.add_argument("-b", "--binsize", default = 0.5, type= float,  help="binsize in m/z units")
parser.add_argument("-w", "--windowsize", help="size of sliding window", default=1, type=float)
parser.add_argument("-p", "--stepsize", help="step size for sliding window", default=0.3, type=float)
parser.add_argument("-m", "--mode", help="mode of operation (1-4)", default=1, type=int)
parser.add_argument("-o", "--output", help="output file name", default = "output.tsv")
args = parser.parse_args()

file = args.fileName
range0 = args.start
range1 = args.end
binsize = args.binsize
numbins = int((range1-range0)/binsize)
window_size = args.windowsize
step_size = args.stepsize
mode = args.mode
output_file = args.output

print("File:    {0:s}".format(args.fileName))
print("Range:   {0:7.1f} to {1:7.1f} m/z".format(args.start, args.end))
print("Binsize: {0:6.1f}".format(args.binsize))
print("Window Size: {0:6.1f} m/z".format(window_size))
print("Step Size:   {0:6.1f} m/z".format(step_size))

pepmass = {}
pepmass = readData(file,range0,range1)

for bin in range(0,numbins):

    lower = float(range0 + bin*binsize)
    upper = float(lower + binsize)

    #print("massbin: {0:12.3f} - {1:12.3f}\t".format(lower,upper), end="")
    count = 0 
    for pep in pepmass.keys():
        mz = float(pepmass[pep])
        if ( mz >= lower and mz < upper):
#            print ("hit: ",pep,mz)
            count += 1
    #print("count:",count)

if mode == 1:
    print("Mode 1: Total count of peptides")
    mode_1_total_count(pepmass, range0, range1, output_file = "mode1.tsv")
elif mode == 2:
    print("Mode 2: Generating histogram")
    mode_2_histogram(pepmass, range0, range1, binsize, numbins, output_file = "mode2.tsv")
elif mode == 3:
    print("Mode 3: Sliding window analysis")
    mode_3_sliding_window(pepmass, range0, range1, window_size, step_size, output_file = "mode3.tsv")
elif mode == 4:
    print("Mode 4: Unique protein identification")
    mode_4_unique_protein(pepmass, output_file = "mode4.tsv")
elif mode == 5:
    print("Mode 5: Fixed_histogram")
    plot_fixed_histogram(file)
elif mode == 6:
    print("Mode 6: Sliding window")
    plot_sliding_window(file)
else:
    print("Invalid mode selected. Please choose between 1 and 4.")

    
if __name__ == "__main__":
    pass
else:
    print("run as module\n")