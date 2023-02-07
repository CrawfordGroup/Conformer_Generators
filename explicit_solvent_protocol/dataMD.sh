#!/bin/bash

# Prompts for user input.
echo -n 'Molecule Name: '
read molecule_name

echo -n 'Spectroscopy: '
read spectroscopy

cwd=$(pwd)

# Defining a function to create and append files.
data () {
    # Finds the number of atoms and assigns it to a variable.
    natom=$(awk '$1 == "NAtoms=" {print $2;exit;}' output.log)
    echo "Number of Atoms: $((natom))"
    number_vibrations=$((3*${natom}-6))
    echo "Number of Vibrations: $((number_vibrations))"
    imaginary_freq=$(awk '$3 == "imaginary" {print $2;exit;}' output.log)
    echo "Number of Imaginary Frequencies: $((imaginary_freq))"
    
    # Finds discrepancy between the number of vibrations and that able to be "read" from a Gaussain output file using awk.
    number_gaussian_freq=$(awk 'BEGIN{count=0;} $1 == "Frequencies" {for(i=3;i<=NF;i++){count+=1}} END{print count}' output.log)
    discrepancy=$((${number_vibrations}-${number_gaussian_freq}))

    # Searches for the word "Frequencies" in the first column and outputs the last three columns into the Freq.txt file.
    #awk '$1 == "Frequencies" { for(i=3;i<=NF;i++){ printf("%.8f\n",$i/=8065.54429)}}' output.log >> Freq.txt
    awk -v imaginary_vibs=$((${imaginary_freq}-${discrepancy})) -v tot_vibs=$((${number_vibrations}-${discrepancy}-1)) '$1 == "Frequencies" {for(i=3;i<=NF;i++){a[j++]=$i}} END{for(j=tot_vibs;j>=imaginary_vibs;j--)printf("%.8f\n",a[j]/=8065.54429)}' output.log >> Freq.txt

    if [ "${spectroscopy}" == "VCD" ]; then
        # Searches for the word "Rot." in the first column and outputs the last three columns into the Intensities.txt file.
        #awk '$1 == "Rot." { for(i=4;i<=NF;i++){ printf("%e\n",$i*=10**-44)}}' output.log >> Intensities.txt
        awk -v imaginary_vibs=${imaginary_freq} -v tot_vibs=$((${number_vibrations}-1)) '$1 == "Rot." {for(i=4;i<=NF;i++){b[j++]=$i}} END{for(j=tot_vibs;j>=imaginary_vibs;j--)printf("%e\n",b[j]*=10**-44)}' output.log >> Intensities.txt
    else
        # Searches for the word "CID3" in the first column and outputs the last three columns into the Intensities.txt file.
        #awk '$1 == "CID3" { for(i=4;i<=NF;i++){ printf("%e\n",$i*=10**4)}}' output.log >> Intensities.txt
        awk -v imaginary_vibs=${imaginary_freq} -v tot_vibs=$((${number_vibrations}-1)) '$1 == "CID3" {for(i=4;i<=NF;i++){b[j++]=$i}} END{for(j=tot_vibs;j>=imaginary_vibs;j--)printf("%e\n",b[j]*=10**4)}' output.log >> Intensities.txt

    fi

    # Prints the snapshot number and free energy into their respective files.
    echo ${conformer_count}  >> $cwd/Conformer.txt
    echo $free_energy >> $cwd/Free_Energy.txt

    # Appends the snapshot number and free energy files together in a columnar manner.
    # Each iteration appends this for calculating the Boltzmann populations later on.
    paste -d'\t' $cwd/Conformer.txt $cwd/Free_Energy.txt >> $cwd/Combined_Free_Energy.txt

    # Removes the excess files.
    rm $cwd/Conformer.txt
    rm $cwd/Free_Energy.txt

    # Creates three new files with the snapshot number and free energy with the same number of lines as will be generated for the frequency and intensities.
    l=1 
    while [ $l -le $((${number_vibrations}-${imaginary_freq})) ];do
        echo ${conformer_count}  >> Conformer.txt
        echo ${free_energy} >> Free_Energy.txt
        let l+=1
    done

    # Appends the snapshot number, free energy, vibrational frequencies, and rotational strengths in a columnar manner.
    paste -d'\t' Conformer.txt Free_Energy.txt Freq.txt Intensities.txt >> Combined_$conformer_count.txt

    # Removes the excess files.
    rm Freq.txt
    rm Intensities.txt
    rm Conformer.txt
    rm Free_Energy.txt

    # Appends the combined file to a master file for further computations.
    cat $cwd/${molecule_name}_MD/cmpd_$conformer_count/Combined_$conformer_count.txt >> $cwd/Combined.txt

    # Removes the excess file.
    rm $cwd/${molecule_name}_MD/cmpd_$conformer_count/Combined_$conformer_count.txt
}                                                                                                                  

echo "Molecule Information"
echo "Molecule Name: ${molecule_name}"
echo ""

CMPD=()
GFE=()
total_conformers=0

number_of_conformers=$(ls -1U "${molecule_name}_MD" | wc -l )
echo "Number of Conformers: $((number_of_conformers))"
echo ""

conformer_count=1
while [ $conformer_count -le $number_of_conformers ]; do           
    echo "Snapshot: $((conformer_count))"

    # Enters the directory for the snapshots.
    cd $cwd/${molecule_name}_MD/cmpd_${conformer_count}

    # Finds the Gibbs free energy and assigns it to a variable.
    free_energy=$(awk -F " " '{if($1=="Sum" && $3=="electronic" && $7=="Energies=") {printf("%.6f",$NF);exit}}' output.log)
    data 

    # Appends the compound and Gibbs free energy arrays.  
    CMPD+=(${conformer_count})
    GFE+=(${free_energy})
    
    echo ""
    let conformer_count+=1
    cd $cwd
done

# Sorts the master file according to vibrational frequency.
sort -nk 3 -o Sorted.txt Combined.txt

# Removes the excess file.
#rm Combined.txt

# Notes:
# First output file "Combined_Free_Energy.txt" is arranged as...
# Snapshot   Gibbs Free Energy (Hartrees)
# Second output file "Sorted.txt" is arranged as...
# Snapshot   Gibbs Free Energy (Hartrees)    Vibrational Frequency (eV)  Rotational Strength (esu^2 cm^2)


