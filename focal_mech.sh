#!/bin/bash

latlow="35.0"
lathig="38.3"
lonlow="21.0"
lonhig="29.0"

rm -rf Focal_mech
mkdir Focal_mech

while read V1 V2 V3 V4 V5 V6 V7 V8 extra 
do
yaer=${V1:6:4}
mnth=${V1:3:2}
daet=${V1:0:2}
howr=${V2:0:2}
minuet=${V2:3:2}
secodn=${V2:6:2}

if [[ "${V7}" -ge 35 ]] && (( $(echo "$V3 >= $latlow" | bc -l) )) && (( $(echo "$V3 <= $lathig" | bc -l) )) && (( $(echo "$V5 >= $lonlow" | bc -l) )) && (( $(echo "$V5 <= $lonhig" | bc -l) )); then
srch=`find /home/user/Documents/Scripts/ -type d -name "${yaer}-${mnth}-${daet}H${howr}M${minuet:0:2}*"`
echo $srch
if [[ "${#srch}" -gt 0 ]]; then
cp -r $srch/ ./Focal_mech/
fi
fi
done < rev_loc.txt

