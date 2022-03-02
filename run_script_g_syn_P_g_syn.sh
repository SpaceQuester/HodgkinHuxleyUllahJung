#!/bin/bash
START=$(date +%s)

for i in $(LC_ALL=C seq 0.0 0.0375 1.50001)
do
	for j in $(LC_ALL=C seq 0.0000 0.0005 0.020001)
	do
        	echo "$i; $j"
        	./UJ_HH_Exp3_LC.out 0.5 $j -0.3 $i 260 0 1
	done
done

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "$(($DIFF / 60)) min $(($DIFF % 60)) sec"
