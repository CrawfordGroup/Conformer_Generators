"""
Created on Mon Feb 27 18:05:06 2023

@author: Aparna K
"""

import os
import pandas as pd
import time

t=time.time()

molecule_name = input("Molecule Name : ")
spectroscopy = input("Spectroscopy : ")

# molecule_name = "A_glucose"
# spectroscopy = "VCD"

cwd=os.getcwd()

print("Molecule Information")
print(f"Molecule Name : {molecule_name}")
print()

CMPD = []
test_list=[]
GFE = []
total_conformers = 0

number_of_conformers = len(os.listdir(f"{molecule_name}_MD"))
print(f"Number of Conformers: {number_of_conformers}")
print()

conformer_count = 1

free_energy_lines = []

while conformer_count <= number_of_conformers:
    freq_lines= []
    intensity_lines = []
    #free_energy_lines = []
    print(f"Snapshot: {conformer_count}")

    # Enters the directory for the snapshots.
    os.chdir(f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}")

    # Finds the Gibbs free energy and assigns it to a variable.
    ctr=0
    with open("output.log", "r") as fout:
        for line in fout:
            splitted_line = line.split()
            ctr+=1
            if len(splitted_line) > 7:
                if splitted_line[0] == "Sum" and splitted_line[2] == "electronic" and splitted_line[6] == "Energies=":
                    #free_energy = float(splitted_line[-1])
                    free_energy = splitted_line[-1]
                    free_energy_lines.append(line)
            if len(splitted_line) > 2:
                if splitted_line[2] == "imaginary":
                    imaginary_freq = int(splitted_line[1])
                if splitted_line[0] == "NAtoms=":
                    natom = int(splitted_line[1])
                if splitted_line[0] == "Frequencies":
                    freq_lines.extend(map(float, splitted_line[2:]))
                if spectroscopy == "VCD":
                    if splitted_line[0] == "Rot.":
                        intensity_lines.extend(map(float, splitted_line[3:]))
                if spectroscopy == "ROA":
                    if splitted_line[0] == "CID3":
                        intensity_lines.extend(map(float, splitted_line[3:]))
        
                
        
                
        print(f"Number of Atoms: {natom}")
        number_vibrations = 3 * natom - 6
        print(f"Number of Vibrations: {number_vibrations}")
        print(f"Number of Imaginary Frequencies: {imaginary_freq}")
    
    # Find discrepancy between the number of vibrations and that able to be "read" from a Gaussian output file using awk.
    
        number_gaussian_freq = len(freq_lines)
    
    #number_gaussian_freq = sum(len(f_line.split()) - 2 for f_line in freq_lines) 
    
    #in bash, we skip the first two fields and starts the counter from field 3, we can avoid a loop by subtracting 2 from the length of line
        discrepancy = number_vibrations - number_gaussian_freq

    # freq_values_we_want = freq_values[tot_vibs:imaginary_vibs-1:-1]
        with open('Freq.txt', 'a') as f:
            imaginary_vibs = imaginary_freq - discrepancy
            tot_vibs = number_vibrations - discrepancy - 1
            for j in range(tot_vibs, imaginary_vibs-1,-1):
                freq = freq_lines[j]
            #for freq in freq_lines[tot_vibs-1:imaginary_vibs-1:-1]:
                f.write(f"{freq/8065.54429:.8f}\n")
                
        with open('Intensities.txt', 'a') as f:
            imaginary_vibs = imaginary_freq
            tot_vibs = number_vibrations - 1
            for j in range(tot_vibs, imaginary_vibs-1,-1):
                intensity = intensity_lines[j]*10**4
                f.write(f"{intensity:.6e}\n")

        
# Writes the conformer count to a file.
    with open(f"{cwd}/Conformer.txt", 'a') as f:
        f.write(f"{conformer_count}\n")
    

# Writes the free energy to a file.
    with open(f"{cwd}/Free_Energy.txt", 'a') as f:
        #f.write(f"{free_energy:.6f}\n")
        #f.write(f"{free_energy:.6f}")
        f.write(free_energy)
        
        
    
# Appends the snapshot number and free energy files together in a columnar manner.
# Each iteration appends this for calculating the Boltzmann populations later on.
        
    conformer_df = pd.read_csv(f"{cwd}/Conformer.txt", sep='\t', header=None)
    free_energy_df = pd.read_csv(f"{cwd}/Free_Energy.txt", sep='\t', header=None)

    combined_df = pd.concat([conformer_df, free_energy_df], axis=1)
    combined_df.to_csv(f"{cwd}/Combined_Free_Energy.txt", sep='\t', index=False, header=False, mode='a')

        
    os.remove(f"{cwd}/Conformer.txt")
    os.remove(f"{cwd}/Free_Energy.txt")

# Creates three new files with the snapshot number and free energy with the same number of lines as will be generated for the frequency and intensities.
    
    with open("Conformer.txt", 'a') as fc, open("Free_Energy.txt", 'a') as fe:
        flag=0
        for l in range(0,number_vibrations-imaginary_freq):
            flag+=1
            fc.write(f"{conformer_count}\n")
            fe.write(f"{free_energy}\n")
    test_list.append(flag)      
    # os.system(cmd)
    conformer_df = pd.read_csv("Conformer.txt", sep='\t', header=None)
    free_energy_df = pd.read_csv("Free_Energy.txt", sep='\t', header=None)
    freq_df = pd.read_csv("Freq.txt", sep='\t', header=None)
    intensities_df = pd.read_csv("Intensities.txt", sep='\t', header=None)

    combined_df = pd.concat([conformer_df, free_energy_df, freq_df, intensities_df], axis=1)
    combined_df.to_csv(f"Combined_{conformer_count}.txt", sep='\t', index=False, header=False, mode='a')
        
    
    # Removes the excess files.
    # os.chdir(f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}")
    # cwd2 = os.getcwd()
    os.remove('Freq.txt')
    os.remove('Intensities.txt')
    os.remove('Conformer.txt')
    os.remove('Free_Energy.txt')
    
    with open(f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/Combined_{conformer_count}.txt", "r") as fcc:
        content = fcc.read()
        with open(f"{cwd}/Combined.txt", "a") as combined_file:
            combined_file.write(content)
                    
    os.remove(f'{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/Combined_{conformer_count}.txt')

    # Appends the compound and Gibbs free energy arrays.
    CMPD.append(conformer_count)
    GFE.append(free_energy)

    print()
    conformer_count += 1
    os.chdir(cwd)

input_file = "Combined.txt"
output_file = "Sorted.txt"


with open(input_file, "r") as finp:
    lines = finp.readlines()

lines = sorted(lines, key=lambda line: float(line.split()[2]))

with open(output_file, "w") as foutp:
    foutp.writelines(lines)
print(time.time()-t)
    
# comb_df = pd.read_csv("Combined", sep='\t', header=None)    
# combined_df.to_csv(f"Combined_{conformer_count}.txt", sep='\t', index=False, header=False, mode='a')






    



