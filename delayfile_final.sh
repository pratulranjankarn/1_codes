#!/bin/bash

rm -f delayfile_final.txt

#awk 'NR > 1' st_delay_husn.txt > test.txt

echo -e "origin_time"'\t' '\t' '\t'"Ev. Lat" '\t'"Ev. Long" '\t'"Ev. Dept" '\t'"Station" "St. Lat" '\t'"St. Long" '\t'"delay 2-4" '\t'"delay 4-8" '\t'"delay 8-16" '\t'"delay 16-32" >> delayfile_final.txt

while read A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11
do
echo $A1
echo -e ${A1}"\t"$A2"\t"$A3"\t"$A4"\t"$A5"\t"$A6"\t"$A7"\t"$A8"\t"$A9"\t"$A10"\t"$A11 >> delayfile_final.txt
done < delayfile_temp.txt
