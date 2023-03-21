import os

molecule_name = str(input("Molecule Name : "))
spectroscopy = str(input("Spectroscopy : "))
cwd = os.getcwd()

compounds = len(os.listdir('{}_MD'.format(molecule_name)))

print("Molecule Name : ",molecule_name)
print("\nSpectroscopy : ",spectroscopy)

cmpd = 1
while cmpd<=compounds:
    os.chdir('{}/{}_MD/cmpd_{}'.format(cwd,molecule_name,cmpd))
    os.system("dos2unix G09_sub_SLURM.sh")
    os.system("sbatch G09_sub_SLURM.sh")
    print("Job submitted for compound ", cmpd)
    cmpd += 1
    os.chdir('{}/{}_MD'.format(cwd,molecule_name))
