#!/bin/bash

latlow="35.0"
lathig="38.3"
lonlow="21.0"
lonhig="29.0"
# run with caution
#: <<'EOF'
#rm -rf Moments
#mkdir Moments
while read V1 V2 V3 V4 V5 V6 V7 V8 extra 
do

#echo $V1 $V2

yaer=${V1:0:4}
mnth=${V1:4:2}
daet=${V1:6:2}
howr=${V2}
rw=${#howr}
if [[ "${rw}" -lt 2 ]]; then
howr="0${howr}"
fi
#howrmin=${V2:3:4}
minuet=${V3}
rw=${#minuet}
if [[ "${rw}" -lt 2 ]]; then
minuet="0${minuet}"
fi
secodn=${V2:6:2}
if [[ "${V7}" -ge 35 ]] && (( $(echo "$V5 >= $latlow" | bc -l) )) && (( $(echo "$V5 <= $lathig" | bc -l) )) && (( $(echo "$V6 >= $lonlow" | bc -l) )) && (( $(echo "$V6 <= $lonhig" | bc -l) )); then
echo ${yaer}-${mnth}-${daet}H${howr}M${minuet}S
fi
qw=`find ./Event_wav_final/ -type d -name "${yaer}-${mnth}-${daet}H${howr}M${minuet}S*"`
echo ${qw}
sw=${#qw}
if [[ "${sw}" -gt 0 ]]; then
wget -qO- http://bbnet.gein.noa.gr/mt_solution/${yaer}/${yaer}${mnth}${daet}_${howr}${minuet}${secodn}.NOA_IG_MT.html | sed -e :a -e 's/<[^>]*>//g;/</N;//ba' > ./Moments/moment_${qw:18:19}.txt

#tail -n +46 new_${V2}.txt > _${V2}.txt

#drnam=$(echo 20${yaer}-${mnth}-${daet}H${howr}M${minuet}S${secodn})
#mkdir $drnam
#rm -f ./${drnam}/Seisgram2k.pick
: << COD
while read M1 M2 M3 M4 M5 M6 M7 M8 M9
do
if [ "${#M1}" -eq "0" ] || [ "${M1}" == "STOP" ]; then

break

fi

if [ "${M4}" != "AML" ] ; then

comp=${M9: -4:1}
if [ "${comp}" = "N" ] ; then
comp="Y"
fi
if [ "${comp}" = "E" ] ; then
comp="X"
fi
#echo $comp
hrpick=$(echo ${M5:0:2})
minpick=$(echo ${M5:3:2})
hrmin=$(echo ${hrpick}${minpick})
secnd=$(echo ${M5:6:6}0)
#echo $secnd
echo $M1 ? $comp ? $M4 ? 20${yaer}${mnth}${daet} ${hrmin} $secnd GAU 0.0 0.0 0.0 0.0 >> ./${drnam}/Seisgram2k.pick
fi
done < pick_${V2}.txt
COD
fi
#rm -f pick_${V2}.txt
#rm -f new_${V2}.txt
done < rev_loc.txt
#EOF
