#! /bin/bash

dirpath="//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/tau_4e-1"
qpath="$dirpath/order_queue.dat"

count=0
while read -r -a line; do

	count=$((count+1))

	if [ $count -ne 1 ]; then

		Ped=${line[0]}
		frameStart=${line[1]}
		frameEnd=${line[2]}
		frameW=${line[3]}
		fname=${line[4]}

		fpathI="$dirpath/Pe$Ped/$fname"
		fpathO="$dirpath/Pe$Ped/lane_order.dat"

		if [ -f $fpathI ]; then

			printf "\n\nProcessing file %d, Ped = %s\n\n" $((count-1)) $Ped

			make order
			./bin/order.o $fpathI $fpathO $frameStart $frameEnd $frameW

		else
			echo "\nInput error! File $fpathI does not exist.\n"
			exit 1
		fi
	fi

done <$qpath

exit 0