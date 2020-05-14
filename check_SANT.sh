#!/bin/bash
# 36.388252
finame="st_delay.txt"

rm -f test.txt test2.txt test3.txt test4.txt

grep "SANT*" $finame > test.txt

awk '{ print $3,$2,$4}' test.txt > test2.txt

while read V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12;
do
echo ${#V9} ${#V10} ${#V11} ${#V12}
grep "Lat $V2" ./loc_Int/husn.sum.grid0.loc.hyp > test3.txt
{
read P1 P2 P3 P4 P5 P6 P7 P8 P9 P10
} < test3.txt
if [[ "${#V9}" -eq 0 ]]; then
V9="0"
fi
if [[ "${#V10}" -eq 0 ]]; then
V10="0"
fi
if [[ "${#V12}" -eq 0 ]]; then
V12="0"
fi
if [[ "${#V11}" -eq 0 ]]; then
V11="0"
fi
echo ${#V9} ${#V10} ${#V11} ${#V12}
grep "Lat $V2" ./loc_Int/egel.sum.grid0.loc.hyp > test3.txt
{
read Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 Q10
} < test3.txt
#echo $Q3/$Q4/${Q5}T${Q6}:${Q7}:${Q8}

if [[ "${#Q1}" -eq 0 ]]; then
dayt=$(echo $P3/$P4/${P5}T${P6}:${P7}:${P8})
echo $P3/$P4/${P5}T${P6}:${P7}:${P8}
secs=`date -d $dayt +%s`
echo $V1 $V2 $V3 $V4 $V5 $V6 $V7 $V8 $V9 $V10 $V11 $V12 $P3 $P4 >> test4.txt
else
dayt=$(echo $Q3/$Q4/${Q5}T${Q6}:${Q7}:${Q8})
echo $Q3/$Q4/${Q5}T${Q6}:${Q7}:${Q8}
secs=`date -d $dayt +%s`
echo $V1 $V2 $V3 $V4 $V5 $V6 $V7 $V8 $V9 $V10 $V11 $V12 $Q3 $Q4 >> test4.txt
fi
done < test.txt
PS1=SANT_plot.ps
REGION="20/29.5/34/39.5"
SCALE=1.1i
GRID="/home/user/Downloads/ETOPO1_Bed_g_gmt4.grd/data"
CPT="custom.cpt"

gmt makecpt -Cgray -T-8000/8550/500 > $CPT
gmt makecpt -Cglobe -T90/150/5 -Z > seis.cpt

rm -f $PS1

#gmt grdgradient $GRID -A225 -Nt0.5 -G$INT

gmt grdimage $GRID -C$CPT -R$REGION -Jm$SCALE -Bf0.5a1/f0.5a1:.:WeSn -V -K -Y1.4 -X2.0 >> $PS1

gmt pscoast -Jm$SCALE -R$REGION -Di -W.3p,0/0/0 -N2p,0/0/0 --DIR_GSHHG=/usr/share/gshhg-gmt-nc4 -V -O -K >> $PS1

gmt psxy test2.txt -R -J -O -V -K -W.1 -Sc.2 -Cseis.cpt -h15 >> $PS1

gmt psscale -D0/3.2/6/1 -B5:Depth:/:km: -Cseis.cpt -O -V >> $PS1

gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tg

rm -f test.txt test2.txt test3.txt
