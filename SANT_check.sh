#!/bin/bash

rm -f test.txt

while read V1 V2 V3 V4 V5 V6 V7 V8 V9;
do
echo $V4 $V2 $V3
grep "Depth $V4*" ../HUSN_EV/loc/husn.20*.hyp | grep "Lat $V2*" | grep "Long $V3*" | awk '{print $3,$4,$5,$6,$7,$8}' >> test.txt
done < SANT_husn_check.txt

paste SANT_husn_check.txt test.txt > SANT_husn_deep.txt
