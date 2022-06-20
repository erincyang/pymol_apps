#!/bin/sh

#$1 should be the pdb to file to read
#split pdb into two files
#the first contains the pdb data
#the second contains the energy data
#also create a res number list from the pdb for renumbering
#the energy resnumbers
		cat $1 | grep '^[A-Z]*_[0-9]\|^[A-Z]*_[A-Z]_[0-9]\|^[A-Z]*:[a-zA-Z]' | rev | cut -d"_" -f1 | rev > energies.dat
		cat $1 | grep "^ATOM " > temp.pdb
		cat temp.pdb | cut -c24-30 | sort -u > res.list
		name=`echo $1 | cut -d"." -f1`

#fix resnumbers
gawk 'FNR==NR{a[NR]=$1;next}{$1=a[FNR]}1' res.list energies.dat > ${1}_energies.txt

mv energies.dat bak.energies.dat
cp ${1}_energies.txt energies.dat

#duplicating rows so that the energies are aligned with the pdb
gawk '
	{
		if(NR == FNR){ #a way of determining that we are reading from the first file
			resi = $1 #residue number
			energy[resi] = $NF #total energy keyed to residue number
		}
		else{ #second file
			e_val=energy[$6] #match key (residue number in pdb) to energy
			print $6, e_val > "energy.dat"
		}
	}
	' energies.dat temp.pdb

#write energy to pdb b-value column
gawk ' 
	NR==FNR { pdb[NR]=$0; next }
	{
	    split(pdb[FNR],flds,FS,seps)
	    flds[11]=$2
	    for (i=1;i in flds;i++)
        	printf "%s%s", flds[i], seps[i]
	    print ""
	}
	' temp.pdb energy.dat > "${name}_energies.pdb"

#cleanup
rm energies.dat
rm energy.dat
rm temp.pdb
rm res.list
rm ${1}_energies.txt
rm bak.energies.dat
