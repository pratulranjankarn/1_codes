#!/bin/bash

dirnam="cr-co_cn_N_orig"
adrsfil="cn_crscr_N_orig"

mkdir $dirnam

while read V1
do
V2=$(echo ${V1:0:55})
filnam=$(echo ${V2: -20:19})
V4=`grep "APE" ${V2}regen.pick | awk '{ print $7,$8,$9}'`
V5=`grep "APE" ${V2}regen.pick | awk '{ print $5}'`
echo ${V5:0:1}
Pyear=${V4:0:4}
Pmnth=${V4:4:2}
Pdate=${V4:6:2}
Pday=`date -d ${Pyear}/${Pmnth}/${Pdate} +%j`
Phr=${V4:9:2}
Pmin=${V4:11:2}
Pcheck=${V4:15:1}
if [[ "$Pcheck" = "." ]]; then
Psec=${V4:(-6):1}
else
Psec=${V4:(-7):2}
fi
Pmsec=${V4:(-4):3}
echo $Pyear
echo $Pday
echo $Phr
echo $Pmin
echo $Psec
echo $Pmsec
echo $Pcheck
sac <<EOF
r $V1
ch A GMT $Pyear $Pday $Phr $Pmin $Psec $Pmsec
wh
cut A 0 30
r
write SAC ./${dirnam}/APE_${filnam}.sacii
write ALPHA ./${dirnam}/APE_${filnam}.ascii
quit
EOF

tail -n +31 ./${dirnam}/APE_${filnam}.ascii > test.txt
awk '{for (i=1;i<=NF;i++) print $i}' test.txt > ./${dirnam}/out_APE_${filnam}.txt

done < ${adrsfil}
