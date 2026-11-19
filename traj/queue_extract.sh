#! /bin/bash

dirpath="//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/tau_1e0_R2"
qpath="$dirpath/coords/particles.dat"

count=0
pnums=""
while read -r -a line; do

	count=$((count+1))

	if [ $count -ne 1 ]; then
		pnums="$pnums${line[0]} "
	fi
		
done <$qpath

fpathI="$dirpath/short_traj.cfg"

if [ -f $fpath ]; then
	g++ -w -g -std=c++17 analysis.cpp extract.cpp -o ./bin/extract.o
	./bin/extract.o $pnums
else
	echo "\nInput error! File $fpathI does not exist.\n"
	exit 1
fi
