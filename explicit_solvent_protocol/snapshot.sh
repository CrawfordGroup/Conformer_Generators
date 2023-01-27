#!/bin/zsh

#___________________
# ==> User Input <==
#___________________

# Prompts for user input.
echo -n 'Molecule Name: '
read molecule_name

echo -n 'Spectroscopy: '
read spectroscopy

echo -n 'Basis Set: '
read basis_set

if [ ${basis_set} = Mixed ]; then
    echo -n 'Solute Basis Set: '
    read solute_basis_set
    echo -n 'Solvent Basis Set: '
    read solvent_basis_set
fi

echo -n 'Snapshots: '
read snapshots

echo -n 'Distance Threshold (Angstroms): '
read D

if [ ${spectroscopy} != VCD -a ${spectroscopy} != ROA ]; then
    echo "${spectroscopy} is not a valid Gaussian spectroscopy supported in this script."
    echo "Valid spectroscopies currently implemented are VCD and ROA."
    exit 1
fi

if [ ${spectroscopy} = ROA ]; then
    echo -n 'Frequency of Incident Radiation (nm): '
    read frequency
fi

cwd=$(pwd)

#_____________________________________________________________
# ==> Generate Data for Single Snapshot from MD Trajectory <==
#_____________________________________________________________

# Read the number of frames in the MD trajectory.
atoms=$(head -1 ${molecule_name}_MD.arc | awk '{print $1}')
lines=$(wc -l ${molecule_name}_MD.arc | awk '{print $1}')
frames=$(($lines/($atoms+2)))
echo "Frames: $frames"
echo " "

# Obtains the specific frames associated with the number of snapshots requested.
line_number=0
counter=1
for snap in {$((${frames}/${snapshots}))..${frames}..$((${frames}/${snapshots}))}; do
    echo "Obtaining coordinate information for frame ${snap}."
    ends=$((${snap}*(${atoms}+2)))
    begins=$((${ends}-${atoms}-1))
    sed -n "${begins},${ends}p;${ends}q" ${molecule_name}_MD.arc > cmpd_${counter}
    let counter+=1
done
echo " "

#__________________________________________________________________________________
# ==> Generate Directories and Modify SLURM Submission and Gaussian Input Files <==
#__________________________________________________________________________________

mkdir "${molecule_name}_MD"
echo "Generating directory structure and modifying input and submission files."

cd "${molecule_name}_MD"

conformer_count=1
while [ $conformer_count -le $snapshots ]; do
    cd $cwd/${molecule_name}_MD/

    # Makes the directory for a specific conformer.
    mkdir "cmpd_$conformer_count"

    # Copies the original submission script into the conformer's directory.
    cp $cwd/G09_sub_SLURM.sh $cwd/${molecule_name}_MD/cmpd_$conformer_count/

    # Modifies the submission script with conformer specific data.
    sed -i '' 's/MOLECULE_NAME_CONFORMER_NUMBER/'${molecule_name}_cmpd_${conformer_count}_${spectroscopy}'/g' $cwd/${molecule_name}_MD/cmpd_$conformer_count/G09_sub_SLURM.sh

    # Makes the submission script executable.
    chmod +x $cwd/${molecule_name}_MD/cmpd_$conformer_count/G09_sub_SLURM.sh

    # Copies the original input file into the conformer's directory.
    cp $cwd/input.dat $cwd/${molecule_name}_MD/cmpd_$conformer_count/

    # Modifies the input file with conformer specific data.
    sed -i '' 's/MOLECULE_NAME_CONFORMER_NUMBER/'${molecule_name}_cmpd_${conformer_count}_${spectroscopy}'/g' $cwd/${molecule_name}_MD/cmpd_$conformer_count/input.dat

    # Modifies the input with spectroscopic specific data.
    sed -i '' 's/SPECTROSCOPY/'${spectroscopy}'/g' $cwd/${molecule_name}_MD/cmpd_$conformer_count/input.dat

    # Modifies the input with basis set specified by the user.
    if [ ${basis_set} != Mixed ]; then
        sed -i '' 's/BASIS_SET/'${basis_set}'/g' $cwd/${molecule_name}_MD/cmpd_$conformer_count/input.dat
    elif [ ${basis_set} = Mixed ]; then
        sed -i '' 's/BASIS_SET/'Gen'/g' $cwd/${molecule_name}_MD/cmpd_$conformer_count/input.dat
    fi

    cd $cwd
    let conformer_count+=1
done
echo " "
cd "${cwd}"

#_____________________________________________________________________
# ==> Modify Snapshots to Include Only Specified Solvent Molecules <==
#_____________________________________________________________________

echo "Generating data for Gaussian input file."

#___________________________
# ==> Initial Parameters <==
#___________________________

# Set solvent atom types.
O=349
H=350

#_____________________________
# ==> Loop Over Compounds  <==
#_____________________________

# Set the file to be read.
for index in {1..${snapshots}}; do
    file=cmpd_${index}
    echo "Compund ${index}"

    #________________________________________________
    # ==> Setup Solute, Solvent, and Total Arrays <==
    #________________________________________________

    # Initializing atom number, atom symbol, X, Y, Z, and atom type arrays for all atoms.
    atom_number=()
    atom_sym=()
    X=()
    Y=()
    Z=()
    atom_type=()

    # Reads the lines in the file and appends to the arrays.
    atom_number+=($(awk '(NR>2){print $1}' cmpd_${index}))
    atom_sym+=($(awk '(NR>2){print $2}' cmpd_${index}))
    X+=($(awk '(NR>2){print $3}' cmpd_${index}))
    Y+=($(awk '(NR>2){print $4}' cmpd_${index}))
    Z+=($(awk '(NR>2){print $5}' cmpd_${index}))
    atom_type+=($(awk '(NR>2){print $6}' cmpd_${index}))

    # Initializes solute X, Y, and Z arrays.
    solute_number=()
    solute_sym=()
    solute_X=()
    solute_Y=()
    solute_Z=()

    # Initializes solvent X, Y, and Z arrays.
    solvent_number=()
    solvent_sym=()
    solvent_X=()
    solvent_Y=()
    solvent_Z=()

    # Appends solute and solvent arrays.
    for i in {1..${#atom_number[@]}}; do
        if [ ${atom_type[i]} -ne $O -a ${atom_type[i]} -ne $H ]; then
            solute_number+=(${atom_number[i]})
            solute_sym+=(${atom_sym[i]})
            solute_X+=(${X[i]})
            solute_Y+=(${Y[i]})
            solute_Z+=(${Z[i]})
        else
            solvent_number+=(${atom_number[i]})
            solvent_sym+=(${atom_sym[i]})
            solvent_X+=(${X[i]})
            solvent_Y+=(${Y[i]})
            solvent_Z+=(${Z[i]})
        fi
    done
    echo "Number of Solute Atoms: ${#solute_number[@]}"

    #___________________________________________
    # ==> Center Solvent Atoms Around Solute <==
    #___________________________________________

    # Read the size of the box in the X, Y, and Z dimensions.

    box_X=$(awk 'NR==2{print $1}' cmpd_${index})
    box_Y=$(awk 'NR==2{print $2}' cmpd_${index})
    box_Z=$(awk 'NR==2{print $3}' cmpd_${index})

    # Compute average coordinate of solute molecule.
    avg_X=0
    avg_Y=0
    avg_Z=0
    for m in {1..${#solute_number[@]}}; do
        let avg_X+=$(bc -l <<< "scale=6; ${solute_X[m]}/${#solute_number[@]}")
        let avg_Y+=$(bc -l <<< "scale=6; ${solute_Y[m]}/${#solute_number[@]}")
        let avg_Z+=$(bc -l <<< "scale=6; ${solute_Z[m]}/${#solute_number[@]}")
    done

    # Set conditions to center the solvent atoms around the solute.
    for a in {1..${#solvent_number[@]}}; do
        if (( $(bc <<< "${solvent_X[a]} - ${avg_X} > ${box_X}/2") )); then
            solvent_X[$a]=$(bc <<< "${solvent_X[$a]} - ${box_X}")
        elif (( $(bc <<< "${solvent_X[a]} - ${avg_X} < -${box_X}/2") )); then
            solvent_X[$a]=$(bc <<< "${solvent_X[$a]} + ${box_X}")
        fi
        if (( $(bc <<< "${solvent_Y[a]} - ${avg_Y} > ${box_Y}/2") )); then
            solvent_Y[$a]=$(bc <<< "${solvent_Y[$a]} - ${box_Y}")
        elif (( $(bc <<< "${solvent_Y[a]} - ${avg_Y} < -${box_Y}/2") )); then
            solvent_Y[$a]=$(bc <<< "${solvent_Y[$a]} + ${box_Y}")
        fi
        if (( $(bc <<< "${solvent_Z[a]} - ${avg_Z} > ${box_Z}/2") )); then
            solvent_Z[$a]=$(bc <<< "${solvent_Z[$a]} - ${box_Z}")
        elif (( $(bc <<< "${solvent_Z[a]} - ${avg_Z} < -${box_Z}/2") )); then
            solvent_Z[$a]=$(bc <<< "${solvent_Z[$a]} + ${box_Z}")
        fi
    done

    #______________________________________________
    # ==> Setup and Compute Intermediate Arrays <==
    #______________________________________________

    # Initialize intermediate solvent arrays for faster computation.
    int_number=()
    int_sym=()
    int_X=()
    int_Y=()
    int_Z=()

    # Compute maximum distance of atoms in solute from average coordinate.
    max_dist=0
    for n in {1..${#solute_number[@]}}; do
        dist=$(bc -l <<< "scale=6; sqrt((${solute_X[n]} - ${avg_X})^2 + (${solute_Y[n]} - ${avg_Y})^2 + (${solute_Z[n]} - ${avg_Z})^2)")
        if (( $(bc -l <<< "$dist>$max_dist") )); then
            max_dist=${dist}
        fi
    done

    # Set the distance threshold for the retaining solvent atoms.
    Dist_thresh=$((${max_dist}+${D}))

    # Compute intermediate solvent arrays.
    for o in {1..${#solvent_number[@]}}; do
        solvent_dist=$(bc -l <<< "scale=6; sqrt((${avg_X} - ${solvent_X[o]})^2 + (${avg_Y} - ${solvent_Y[o]})^2 + (${avg_Z} - ${solvent_Z[o]})^2)")
        if (( $(bc -l <<< "${solvent_dist}<=${Dist_thresh}") )); then
            int_number+=(${solvent_number[o]})
            int_sym+=(${solvent_sym[o]})
            int_X+=(${solvent_X[o]})
            int_Y+=(${solvent_Y[o]})
            int_Z+=(${solvent_Z[o]})
        fi
    done
    echo "Number of Intermediate Solvent Atoms: ${#int_number[@]}"

    #________________________________________________
    # ==> Setup and Compute Gaussian Input Arrays <==
    #________________________________________________

    # Initialize arrays for Gaussian input.
    final_number=()
    final_sym=()
    final_X=()
    final_Y=()
    final_Z=()

    # Appending arrays for Gaussian input with solute data.
    final_number+=(${solute_number[@]})
    final_sym+=(${solute_sym[@]})
    final_X+=(${solute_X[@]})
    final_Y+=(${solute_Y[@]})
    final_Z+=(${solute_Z[@]})

    # Determines solvent molecules distance from solute atoms and appends.
    for j in {1..${#solute_number[@]}}; do
        #if [ ${solute_sym[j]} != H ]; then
            for k in {1..${#int_number[@]}}; do
                d=$(bc -l <<< "scale=6; sqrt((${solute_X[j]} - ${int_X[k]})^2 + (${solute_Y[j]} - ${int_Y[k]})^2 + (${solute_Z[j]} - ${int_Z[k]})^2)")
                if (( $(bc -l <<< "$d<=$D") )); then
                    no_duplicates=true
                    for l in {1..${#final_number[@]}}; do
                        if [ ${int_number[k]} -eq ${final_number[l]} ]; then
                            no_duplicates=false
                        fi
                    done
                    if ${no_duplicates}; then
                        final_number+=(${int_number[k]})
                        final_sym+=(${int_sym[k]})
                        final_X+=(${int_X[k]})
                        final_Y+=(${int_Y[k]})
                        final_Z+=(${int_Z[k]})
                    fi
                fi
            done
        #fi
    done
    echo "Number of Atoms before Solvent Check: ${#final_number[@]}"

    #___________________________________________________
    # ==> Confirm Presence of Full Solvent Molecules <==
    #___________________________________________________

    # Determine if all solvent atoms are part of a full molecule.
    solvent_complete=false
    while [ ${solvent_complete} = false ]; do

        # Initializing solvent connectivity arrays.
        connect_1=()
        connect_2=()
        connect_3=()
        connect_4=()
        col=()

        # Read solvent specific lines of the file and append to the arrays.
        for p in {1..${#final_number[@]}}; do
            ln_num=${final_number[p]}
            columns=$(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print NF}' cmpd_${index})
            if [ ${columns} -eq 6 ]; then
                connect_1+=(0)
                connect_2+=(0)
                connect_3+=(0)
                connect_4+=(0)
            elif [ ${columns} -eq 7 ]; then
                connect_1+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $7}' cmpd_${index}))
                connect_2+=(0)
                connect_3+=(0)
                connect_4+=(0)
            elif [ ${columns} -eq 8 ]; then
                connect_1+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $7}' cmpd_${index}))
                connect_2+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $8}' cmpd_${index}))
                connect_3+=(0)
                connect_4+=(0)
            elif [ ${columns} -eq 9 ]; then
                connect_1+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $7}' cmpd_${index}))
                connect_2+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $8}' cmpd_${index}))
                connect_3+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $9}' cmpd_${index}))
                connect_4+=(0)
            elif [ ${columns} -eq 10 ]; then
                connect_1+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $7}' cmpd_${index}))
                connect_2+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $8}' cmpd_${index}))
                connect_3+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $9}' cmpd_${index}))
                connect_4+=($(awk -v ln="${ln_num}" '$1==ln && $2!="Great" {print $10}' cmpd_${index}))
            fi
            col+=(${columns})
        done

        # Initialize new atom array.
        new_number=()

        # Checking status of atoms connectivity.
        for q in {1..${#final_number[@]}}; do

            # Initializing status of connectivity.
            connect_1_present=false
            connect_2_present=false
            connect_3_present=false
            connect_4_present=false
            for r in {1..${#connect_1[@]}}; do
                if [ ${connect_1[q]} -eq ${final_number[r]} ]; then
                    connect_1_present=true
                elif [ ${connect_2[q]} -eq ${final_number[r]} ]; then
                    connect_2_present=true
                elif [ ${connect_3[q]} -eq ${final_number[r]} ]; then
                    connect_3_present=true
                elif [ ${connect_4[q]} -eq ${final_number[r]} ]; then
                    connect_4_present=true
                fi
                if [ ${connect_1[q]} -eq 0 ]; then
                    connect_1_present=true
                    connect_2_present=true
                    connect_3_present=true
                    connect_4_present=true
                elif [ ${connect_2[q]} -eq 0 ]; then
                    connect_2_present=true
                    connect_3_present=true
                    connect_4_present=true
                elif [ ${connect_3[q]} -eq 0 ]; then
                    connect_3_present=true
                    connect_4_present=true
                elif [ ${connect_4[q]} -eq 0 ]; then
                    connect_4_present=true
                fi
            done

            # Checking for connectivity and duplicates.
            if [[ ${connect_1_present} = false ]]; then
                no_duplicates=true
                for t in {1..${#new_number[@]}}; do
                    if [[ ${connect_1[q]} -eq ${new_number[t]} ]]; then
                        no_duplicates=false
                    fi
                done
                if ${no_duplicates}; then
                    new_number+=(${connect_1[q]})
                fi
            elif [[ ${connect_2_present} = false ]]; then
                no_duplicates=true
                for t in {1..${#new_number[@]}}; do
                    if [[ ${connect_2[q]} -eq ${new_number[t]} ]]; then
                        no_duplicates=false
                    fi  
                done
                if ${no_duplicates}; then
                    new_number+=(${connect_2[q]})
                fi  
            elif [[ ${connect_3_present} = false ]]; then
                no_duplicates=true
                for t in {1..${#new_number[@]}}; do
                    if [[ ${connect_3[q]} -eq ${new_number[t]} ]]; then
                        no_duplicates=false
                    fi  
                done
                if ${no_duplicates}; then
                    new_number+=(${connect_3[q]})
                fi  
            elif [[ ${connect_4_present} = false ]]; then
                no_duplicates=true
                for t in {1..${#new_number[@]}}; do
                    if [[ ${connect_4[q]} -eq ${new_number[t]} ]]; then
                        no_duplicates=false
                    fi  
                done
                if ${no_duplicates}; then
                    new_number+=(${connect_4[q]})
                fi  
            fi
        done
        echo "Number of New Atoms: ${#new_number[@]}"
        if [ ${#new_number[@]} -eq 0 ]; then
            echo "Solvent Complete."
            solvent_complete=true
            break
        else
            echo "Solvent Incomplete. Appending and obtaining connectivity of new atoms."
        fi

        # Appending new atoms to final atom arrays.
        for s in {1..${#new_number[@]}}; do
            ind=$((${new_number[s]}-${#solute_number[@]}))
            final_number+=(${solvent_number[${ind}]})
            final_sym+=(${solvent_sym[${ind}]})
            final_X+=(${solvent_X[${ind}]})
            final_Y+=(${solvent_Y[${ind}]})
            final_Z+=(${solvent_Z[${ind}]})
        done
    done
    echo "Final Number of Atoms: ${#final_number[@]}"
    echo " "

    #______________________________________________
    # ==> Writing Guassian Input Arrays to File <==
    #______________________________________________

    # Overwriting compound files with final data.
    for i in {1..${#final_number[@]}}; do
        printf "%-2s \t %-10s \t %-10s \t %-10s\n" "${final_sym[i]}" "${final_X[i]}" "${final_Y[i]}" "${final_Z[i]}" >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
    done

    echo " " >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
    
    if [ ${basis_set} = Mixed ]; then
        echo "1 - ${#solute_number[@]} 0" >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
        echo "${solute_basis_set}" >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
        echo "****" >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
        echo "$((${#solute_number[@]}+1)) - ${#final_number[@]} 0" >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
        echo "${solvent_basis_set}" >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
        echo "****" >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
    fi  

    echo " " >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
    echo " " >> $cwd/${molecule_name}_MD/cmpd_${index}/input.dat
done

cd ${cwd}
rm cmpd_*


