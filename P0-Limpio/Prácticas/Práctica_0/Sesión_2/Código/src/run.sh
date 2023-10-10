#!/bin/bash

make sec

in_file="../input/in.txt"
out_file="../output/out.txt"

while read line ; do 

	if [[ $line != \#* ]]
	then
		n=$(echo $line | tr -s ' ' | cut -f1 -d ' ')
		gen=$(echo $line | tr -s ' ' | cut -f2 -d ' ')
		tam=$(echo $line | tr -s ' ' | cut -f3 -d ' ')
	
		echo -e >> $out_file
		echo -n "Executing with: " >> $out_file
		echo -e "N = "$n" N_GEN = "$gen" TAM_POB = "$tam >> $out_file
		make test_sec N=$n N_GEN=$gen T_POB=$tam
	fi
done < $in_file
