# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 09:58:32 2023

@author: Aparna K
"""

import numpy as np
import h5py
import time 
import os
import shutil
import math

# User Input prompts
t3=time.time()
molecule_name = input('Molecule Name: ')
spectroscopy = input('Spectroscopy: ')
basis_set = input('Basis Set: ')

if basis_set == 'Mixed':
    solute_basis_set = input('Solute Basis Set: ')
    solvent_basis_set = input('Solvent Basis Set: ')

snapshots = input('Number of Snapshots: ')
snapshots = int(snapshots)
D = input('Distance Threshold (Angstroms): ')

if spectroscopy != 'VCD' and spectroscopy != 'ROA':
    print(f"{spectroscopy} is not a valid Gaussian spectroscopy supported in this script.")
    print("Valid spectroscopies currently implemented are VCD and ROA.")
    exit(1)

if spectroscopy == 'ROA':
    frequency = input('Frequency of Incident Radiation (nm): ')

cwd = os.getcwd()

indx =[]
coord = []
at_type = []
connect = []
dims_box = np.zeros((3))
frame_no = 1

t1 = time.time()
with h5py.File(f"{molecule_name}_MD.h5" , "w") as h5:
    with open(f"{molecule_name}_MD.arc" , "r") as f:
        content = f.readlines()
        lines = len(content)
        content.pop(0)
        for i , line in enumerate(content):
            splt = line.split()
            if splt[1] != "Great":
                if splt[-1] == "90.000000":
                    for a in range(3):
                        dims_box[a] = float(splt[a])
                elif len(splt)>6:
                    connect_line = np.zeros((4))
                    xyz_line = np.zeros((3))
                    indx.append(int(splt[0]))
                    at_type.append(int(splt[5]))

                    for a in range(3):
                        xyz_line[a]  = float(splt[a+2])
                    coord.append(xyz_line)

                    for a in range(len(splt)-6):
                        connect_line[a] = int(splt[6+a])
                    connect.append(connect_line)
            
            elif splt[1]=="Great" or not line:
                grp = h5.create_group(f"frame{frame_no}")
                grp.create_dataset("dims_box" , data = np.array(dims_box))
                grp.create_dataset("indx" ,     data = np.array(indx))
                grp.create_dataset("coord" ,    data = np.array(coord))
                grp.create_dataset("at_type" ,  data = np.array(at_type))
                grp.create_dataset("connect" ,  data = np.array(connect))
                frame_no +=1
                no_atoms = []
                indx =[]
                at_sym = []
                coord = []
                at_type = []
                connect = []
                dims_box = np.zeros((3))
        grp = h5.create_group(f"frame{frame_no}")
        grp.create_dataset("dims_box" , data = np.array(dims_box))
        grp.create_dataset("indx" ,     data = np.array(indx))
        grp.create_dataset("coord" ,    data = np.array(coord))
        grp.create_dataset("at_type" ,  data = np.array(at_type))
        grp.create_dataset("connect" ,  data = np.array(connect))

t2 = time.time()
print(t2-t1)

with open(f"{molecule_name}_MD.arc","r") as f:
      atoms = int(f.readline().split()[0])
    
frames = (lines // (atoms + 2))  #2 extra to account for GROMACS separator and box dimensions
print(f"Frames: {frames}")
print()


#Generating Directories and Modifying SLURM submission and Gaussian Input files
os.mkdir(f"{molecule_name}_MD")
print("Generating directory structure and modifying input and submission files.")
os.chdir(f"{molecule_name}_MD")

conformer_count = 1
while conformer_count <= snapshots:

    # Makes the directory for a specific conformer.
    os.mkdir(f"cmpd_{conformer_count}")

    # Copies the original submission script into the conformer's directory.
    shutil.copy(f"{cwd}/G09_sub_SLURM.sh", f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/")

    # Makes the submission script executable.
    os.chmod(f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/G09_sub_SLURM.sh", 0o755)

    # Copies the original input file into the conformer's directory.
    shutil.copy(f"{cwd}/input.dat", f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/")

    # Modifies the input file with conformer, spectroscopy specific data and inputs the basis set specified by the user.
    with open(f"{cwd}/{molecule_name}_MD/cmpd_{conformer_count}/input.dat", "r+") as file:
        content = file.read()
        content = content.replace("MOLECULE_NAME_CONFORMER_NUMBER", f"{molecule_name}_cmpd_{conformer_count}_{spectroscopy}")
        content = content.replace("SPECTROSCOPY", spectroscopy)
        content = content.replace("BASIS_SET", basis_set if basis_set != "Mixed" else "Gen")
        file.seek(0)
        file.write(content)
        file.truncate()

    conformer_count += 1

print()
os.chdir(cwd)

#Modifying Snapshots to include only specified solvent molecules

print("Generating data for Gaussian input file.")

#Initial params

# Set solvent atom types.
O=349
H=350

flag = frames//int(snapshots)
temp_indx = flag

# Set the file to be read.
for index in range(1, snapshots+1):
#with h5py.File(f"{molecule_name}_MD.h5", "r") as f:
    #flag = frames//int(snapshots)
    #temp_indx = flag
    with h5py.File(f"{molecule_name}_MD.h5", "r") as f:
    #for index in range(1, snapshots+1):
        print(f"Obtaining coordinate information for frame {temp_indx} \n")
        atom_sym = []
        intmdt = f"frame{temp_indx}/"
        box = np.array(f[intmdt]['dims_box'])
        atom_number = np.array(f[intmdt]['indx'])
        coords = np.array(f[intmdt]['coord'])
        atom_type = np.array(f[intmdt]['at_type'])
        connect_array = np.array(f[intmdt]['connect'])
        for i in range(len(connect_array)):
            if connect_array[i][0] and connect_array[i][1] and connect_array[i][2] and connect_array[i][3]:
                atom_sym.append('C')
            elif connect_array[i][0] and connect_array[i][1]:
                atom_sym.append('O')
            elif connect_array[i][0]:
                atom_sym.append('H')
    
        
        # Initializes solute X, Y, and Z arrays.
        solute_number = []
        solute_sym = []
        solute_X = []
        solute_Y = []
        solute_Z = []
        
        # Initializes solvent X, Y, and Z arrays.
        solvent_number = []
        solvent_sym = []
        solvent_X = []
        solvent_Y = []
        solvent_Z = []
        
        # Appends solute and solvent arrays.
        for i in range(len(atom_number)):
            if atom_type[i] != O and atom_type[i] != H:
                solute_number.append(atom_number[i])
                solute_sym.append(atom_sym[i])
                solute_X.append(coords[i][0])
                solute_Y.append(coords[i][1])
                solute_Z.append(coords[i][2])
            else:
                solvent_number.append(atom_number[i])
                solvent_sym.append(atom_sym[i])
                solvent_X.append(coords[i][0])
                solvent_Y.append(coords[i][1])
                solvent_Z.append(coords[i][2])
        print(f'Number of Solute Atoms: {len(solute_number)}')
        
        # Compute average coordinate of solute molecule.
        avg_X = np.mean(solute_X)
        avg_Y = np.mean(solute_Y)
        avg_Z = np.mean(solute_Z)
        
        # Set conditions to center the solvent atoms around the solute.
        for a in range(len(solvent_number)):
            if (solvent_X[a]-avg_X) > box[0]/2:
                solvent_X[a] -= box[0]
            elif (solvent_X[a]-avg_X) < -box[0]/2:
                solvent_X[a] += box[0]
            if (solvent_Y[a]-avg_Y) > box[1]/2:
                solvent_Y[a] -= box[1]
            elif (solvent_Y[a]-avg_Y) < -box[1]/2:
                solvent_Y[a] += box[1]
            if (solvent_Z[a]-avg_Z) > box[2]/2:
                solvent_Z[a] -= box[2]
            elif (solvent_Z[a]-avg_Z) < -box[2]/2:
                solvent_Z[a] += box[2]
         
        
        #Computing inetermediate arrays
        
        # Initialize intermediate solvent arrays for faster computation.
        int_number=[]
        int_sym=[]
        int_X=[]
        int_Y=[]
        int_Z=[]
        
        # Compute maximum distance of atoms in solute from average coordinate.
        max_dist = 0
        for n in range(len(solute_number)):
            dist = math.sqrt((solute_X[n] - avg_X)**2 + (solute_Y[n] - avg_Y)**2 + (solute_Z[n] - avg_Z)**2)
            if dist > max_dist:
                max_dist = dist
        
        # Set the distance threshold for the retaining solvent atoms.
        Dist_thresh = max_dist + float(D)
        
        # Compute intermediate solvent arrays.
        for o in range(len(solvent_number)):
            solvent_dist = math.sqrt((avg_X - solvent_X[o])**2 + (avg_Y - solvent_Y[o])**2 + (avg_Z - solvent_Z[o])**2)
            if solvent_dist <= Dist_thresh:
                int_number.append(solvent_number[o])
                int_sym.append(solvent_sym[o])
                int_X.append(solvent_X[o])
                int_Y.append(solvent_Y[o])
                int_Z.append(solvent_Z[o])
            
        
        print(f"Number of Intermediate Solvent Atoms: {len(int_number)}")
        
        # ==> Setup and Compute Gaussian Input Arrays <==
        
        # Initialize arrays for Gaussian input.
        final_number = []
        final_sym = []
        final_X = []
        final_Y = []
        final_Z = []
        
        # Appending arrays for Gaussian input with solute data.
        final_number.extend(solute_number)
        final_sym.extend(solute_sym)
        final_X.extend(solute_X)
        final_Y.extend(solute_Y)
        final_Z.extend(solute_Z)
        
        # Determines solvent molecules distance from solute atoms and appends.
        for j in range(len(solute_number)):
            for k in range(len(int_number)):
                d = math.sqrt((solute_X[j] - int_X[k])**2 + (solute_Y[j] - int_Y[k])**2 + (solute_Z[j] - int_Z[k])**2)
                if d <= float(D):
                    no_duplicates = True
                    for l in range(len(final_number)):
                        if int_number[k] == final_number[l]:
                            no_duplicates = False
                            break #doesn't matter
                    if no_duplicates:
                        final_number.append(int_number[k])
                        final_sym.append(int_sym[k])
                        final_X.append(int_X[k])
                        final_Y.append(int_Y[k])
                        final_Z.append(int_Z[k])
        
        print(f"Number of Atoms before Solvent Check: {len(final_number)}")
        
        #Confirming the presence of full solvent molecules
        
        solvent_complete = False
        poorni = 0
        
        while not solvent_complete:
                    
            # # Initialize new atom array.
            new_number = []
            con = []
            mod_connect = np.column_stack((atom_number, connect_array))
            temp_ctr = 0
            for p in range(len(final_number)):
                ln_num = final_number[p]
                for q in range(len(atom_number)):
                    if ln_num == atom_number[q]:
                        for g in range(len(connect_array)):
                            if mod_connect[g][0] == ln_num:
                                con.append(mod_connect[g][1:]) 
                                break
            con = np.array(con)
            connect = [con[:,0], con[:,1], con[:,2], con[:,3]]
            
            #Checking status of atoms connectivity
            for q in range(len(final_number)):
                lb= [False, False, False, False]
                for r in range(len(connect[0])):
                    for y in range(4):
                        if connect[y][q] == final_number[r]:
                            lb[y] = True
                            break
                    for y in range(4):
                        if connect[y][q] == 0:
                            if y<4:
                                lb[y] = True
                            if y<3:
                                lb[y] = True
                            if y<2:
                                lb[y] = True
                            if y<1:
                                lb[y] = True
                for y in range(4):
                    if lb[y] == False:
                        no_duplicates = True
                        for t in range(len(new_number)):
                            #print(len(new_number))
                            if connect[y][q] == new_number[t]:
                                no_duplicates = False
                                break
                        if no_duplicates:
                            new_number.append(connect[y][q])
                            break
                    
            
                
            # Checking if solvent is complete.
            print("Number of New Atoms:", len(new_number))
            if not len(new_number):
                print("Solvent Complete.")
                solvent_complete = True
                break
            else:
                print("Solvent Incomplete. Appending and obtaining connectivity of new atoms.")
        
            # Appending new atoms to final atom arrays.
            for s in range(len(new_number)):
                ind = new_number[s] - len(solute_number) -1
                ind = int(ind)
                final_number.append(solvent_number[ind])
                final_sym.append(solvent_sym[ind])
                final_X.append(solvent_X[ind])
                final_Y.append(solvent_Y[ind])
                final_Z.append(solvent_Z[ind])
            poorni += 1
            
        print("Final Number of Atoms:", len(final_number))
        print(" ")
        
        #Writing Gaussian Input Arrays to file
          
        # Overwriting compound files with final data.
        for i in range(len(final_number)):
            with open(os.path.join(cwd, f"{molecule_name}_MD", f"cmpd_{index}", "input.dat"), "a") as f:
                f.write("{:<2} \t {:<10} \t {:<10} \t {:<10}\n".format(final_sym[i], final_X[i], final_Y[i], final_Z[i]))
        
        with open(os.path.join(cwd, f"{molecule_name}_MD", f"cmpd_{index}", "input.dat"), "a") as f:
            f.write("\n")   
        
            if basis_set == "Mixed":
                with open(os.path.join(cwd, f"{molecule_name}_MD", f"cmpd_{index}", "input.dat"), "a") as f:
                    f.write("1 - {} 0\n".format(len(solute_number)))
                    f.write("{}\n".format(solute_basis_set))
                    f.write("****\n")
                    f.write("{} - {} 0\n".format(len(solute_number)+1, len(final_number)))
                    f.write("{}\n".format(solvent_basis_set))
                    f.write("****\n")
                    
        with open(os.path.join(cwd, f"{molecule_name}_MD", f"cmpd_{index}", "input.dat"), "a") as f:
            f.write(f"{frequency} nm")
            f.write("\n")
            f.write("\n")
            f.write("Surface=SAS")
            f.write("\n")
            f.write("\n")
            
    temp_indx += flag
            
    

            
        
os.chdir(cwd)
#no need to remove intermediate files in hdf5 formalism!!
# for index in range(1,snapshots+1):
#     os.remove(f"cmpd_{index}")
t2=time.time()
print(t2-t3)




mod_connect = np.column_stack((atom_number, connect_array))
idx = np.isin(atom_number, final_number)
con = mod_connect[idx][:, 1:]
connect = np.transpose(con)



