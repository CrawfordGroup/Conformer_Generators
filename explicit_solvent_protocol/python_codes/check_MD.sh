#!/bin/bash

echo -n "Molecule Name: "
read molecule_name

echo -n "Spectroscopy: "
read spectroscopy

cwd=$(pwd)

cd $cwd

echo "Molecule Name: ${molecule_name}"

number_of_conformers=$(ls -1U "${molecule_name}_MD" | wc -l )
echo "Number of Conformers: $((number_of_conformers))"

conformer_count=1
while [ $conformer_count -le $number_of_conformers ]; do
        cd $cwd/${molecule_name}_MD/cmpd_${conformer_count}
        echo "Conformer: $((conformer_count))"
        if [ `grep -cr "Normal termination of Gaussian 09" $cwd/${molecule_name}_MD/cmpd_${conformer_count}/output.log` = 1 ]; then
                echo "Status: Complete"
        else
                echo "Status: Incomplete"
        fi
        let conformer_count+=1
done
