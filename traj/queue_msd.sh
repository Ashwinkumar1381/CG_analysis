#! /bin/bash

dirpath="//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100"
qpath="$dirpath/short_msd_queue.dat"

count=0
while read -r -a line; do

	count=$((count+1))

	if [ $count -ne 1 ]; then

		Ped=${line[0]}
		frameStart=${line[1]}
		frameEnd=${line[2]}
		frameW=${line[3]}
		fname=${line[4]}

		fpathI="$dirpath/tau_$Ped/$fname"
		fpathO="$dirpath/tau_$Ped/short_msd.dat"

		if [ -f $fpathI ]; then

			printf "\n\nProcessing file %d, Ped = %s\n\n" $((count-1)) $Ped

			make msd
			./bin/msd.o $fpathI $fpathO $frameStart $frameEnd $frameW

		else
			echo "\nInput error! File $fpathI does not exist.\n"
			exit 1
		fi
	fi

done <$qpath

exit 0