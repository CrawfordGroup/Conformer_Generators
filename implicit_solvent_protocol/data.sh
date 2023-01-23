#!/bin/bash

# Prompts for user input.
echo -n 'Molecule Name: '
read molecule_name

echo -n 'Spectroscopy: '
read spectroscopy

ring_structures=$(ls -1U "${molecule_name}" | wc -l )
cwd=$(pwd)

# Defining a function to create and append files.
data () {
    # Searches for the word "Frequencies" in the first column and outputs the last three columns into the Freq.txt file.
    awk '$1 == "Frequencies" { for(i=3;i<=NF;i++){ printf("%.8f\n",$i/=8065.54429)}}' output.log >> Freq.txt
    
    if [ "${spectroscopy}" == "VCD" ]; then
        # Searches for the word "Rot." in the first column and outputs the last three columns into the RotStr.txt file.
        awk '$1 == "Rot." { for(i=4;i<=NF;i++){ printf("%e\n",$i*=10**-44)}}' output.log >> Intensities.txt
    else
        # Searches for the word "CID3" in the first column and outputs the last three columns into the RotStr.txt file.
        awk '$1 == "CID3" { for(i=4;i<=NF;i++){ printf("%e\n",$i*=10**4)}}' output.log >> Intensities.txt
    fi
    # Finds the number of atoms and assigns it to a variable.
    natom=$(awk '$1 == "NAtoms=" {print $2;exit;}' output.log)

    # Prints the ring structure, conformer number, and free energy into their respective files.
    echo ${rs} >> $cwd/Ring_Structure.txt
    echo ${conformer_count}  >> $cwd/Conformer.txt
    echo $free_energy >> $cwd/Free_Energy.txt

    # Appends the ring structure, conformer number, and free energy files together in a columnar manner.
    # Each iteration appends this for calculating the Boltzmann populations later on.
    paste -d'\t' $cwd/Ring_Structure.txt $cwd/Conformer.txt $cwd/Free_Energy.txt >> $cwd/Combined_Free_Energy.txt

    # Removes the excess files.
    rm $cwd/Ring_Structure.txt
    rm $cwd/Conformer.txt
    rm $cwd/Free_Energy.txt

    # Creates three new files with the ring structure, conformer count, and free energy with the same number of lines as will be generated for the frequency and intensities.
    number_vibrations=$((3*${natom}-6))
    l=1 
    while [ $l -le ${number_vibrations} ];do
        echo ${rs} >> Ring_Structure.txt
        echo ${conformer_count}  >> Conformer.txt
        echo ${free_energy} >> Free_Energy.txt
        let l+=1
    done

    # Appends the ring structure, conformer count, free energy, vibrational frequencies, and rotational strengths in a columnar manner.
    paste -d'\t' Ring_Structure.txt Conformer.txt Free_Energy.txt Freq.txt Intensities.txt >> Combined_$conformer_count.txt

    # Removes the excess files.
    rm Freq.txt
    rm Intensities.txt
    rm Conformer.txt
    rm Free_Energy.txt
    rm Ring_Structure.txt

    # Appends the combined file to a master file for further computations.
    cat $cwd/${molecule_name}/rs_${rs}/cmpd_$conformer_count/Combined_$conformer_count.txt >> $cwd/Combined.txt

    # Removes the excess file.
    rm $cwd/${molecule_name}/rs_${rs}/cmpd_$conformer_count/Combined_$conformer_count.txt
}                                                                                                                  

echo "Molecule Information"
echo "Molecule Name: ${molecule_name}"
echo "Number of Ring Structures: $((ring_structures))"
echo ""

RS=()
CMPD=()
GFE=()
total_conformers=0
rs=1

while [ $rs -le $ring_structures ]; do
    number_of_conformers=$(ls -1U "${molecule_name}/rs_${rs}" | wc -l )
    echo "Ring Structure: $((rs))"
    echo "Number of Conformers: $((number_of_conformers))"
    echo ""

    conformer_count=1
    while [ $conformer_count -le $number_of_conformers ]; do           
        echo "Ring Structure: $((rs))"
        echo "Conformer: $((conformer_count))"
        # Enters the directory for the ring structure and conformer.
        cd $cwd/${molecule_name}/rs_${rs}/cmpd_${conformer_count}

        # Finds the Gibbs free energy and assigns it to a variable.
        free_energy=$(awk -F " " '{if($1=="Sum" && $3=="electronic" && $7=="Energies=") {printf("%.6f",$NF);exit}}' output.log)
        Unique=true
        if [ "${#RS[@]}" == 0 ]; then
            echo "Unique: ${Unique}"
            data
            RS+=(${rs})
            CMPD+=(${conformer_count})
            GFE+=(${free_energy})
        else
            i=0
            while [ $i -lt ${#GFE[@]} ]; do
                if (( $(echo "${free_energy} != ${GFE[i]}" |bc -l) )) ; then
                    #echo -e $i '\t' ${free_energy} '\t' ${GFE[i]}
                    let i+=1
                else
                    Unique=false
                    let i=${#GFE[@]}
                fi
            done
                if [ "${Unique}" = true ]; then
                    echo "Unique: ${Unique}"
                    data
                    # Appends the ring structure, compound, and Gibbs free energy arrays.  
                    RS+=(${rs})
                    CMPD+=(${conformer_count})
                    GFE+=(${free_energy})
                else
                    echo "Unique: ${Unique}"
                fi
        fi
        echo ""
        let conformer_count+=1
        cd $cwd
    done
    let rs+=1
done

# Sorts the master file according to vibrational frequency.
sort -nk 4 -o Sorted.txt Combined.txt

# Removes the excess file.
rm Combined.txt

# Notes:
# First output file "Combined_Free_Energy.txt" is arranged as...
# Ring Structure    Conformer   Gibbs Free Energy (Hartrees)
# Second output file "Sorted.txt" is arranged as...
# Ring Structure    Conformer   Gibbs Free Energy (Hartrees)    Vibrational Frequency (eV)  Rotational Strength (esu^2 cm^2)


