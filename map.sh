#!/bin/bash
: <<SOD
rm -f del24.txt del48.txt del816.txt del1632.txt

pi=$(echo "scale=10; 4*a(1)" | bc -l)

function asin() {
    if (( $(echo "$1 == 1" | bc -l) ));then
       echo "90"   
    elif (( $(echo "$1 < 1" | bc -l) ));then
       echo "scale=15;a(sqrt((1/(1-($1^2)))-1))" | bc -l
    elif (( $(echo "$1 > 1" | bc -l) ));then
       echo "error"
    fi
}

awk 'NR>=2' st_delay.txt > test.txt
while read V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12
do
echo $V1 $V4
if (( $(echo "${V4} > 40" | bc -l) )); then
olatrad=$(echo "${V2}*${pi}/180" | bc -l)
olonrad=$(echo "${V3}*${pi}/180" | bc -l)
slatrad=$(echo "${V5}*${pi}/180" | bc -l)
slonrad=$(echo "${V6}*${pi}/180" | bc -l)

#echo $olatrad $olonrad $slatrad $slonrad

delon=$(echo "${olonrad}-${slonrad}" | bc)
delat=$(echo "${olatrad}-${slatrad}" | bc)

#echo $delon $delat

alpha_1=$(echo "s(${delat}/2)^2" | bc -l)
alpha_2=$(echo "c(${olatrad})*c(${slatrad})*(s(${delon}/2)^2)" | bc -l)

#echo ${alpha_1} ${alpha_2}

alpha=$(echo "sqrt(${alpha_1}+${alpha_2})" | bc -l)
beta=`asin ${alpha}`

#echo $alpha $beta

stdist=$(echo "scale=5;6371*2*${beta}" | bc -l)

dist=$(echo "scale=5;sqrt(${stdist}^2+${V4}^2)" | bc -l)
#echo $dist $V4
if (( $(echo "50 <= ${dist} && ${dist} <= 250" | bc -l) )); then
#echo $dist
if (( $(echo "0.6 < ${V7} && ${V7} <= 10" | bc -l) )) && [[ "${#V7}" -gt 0 ]]; then
#echo $V6 $V5 $V7
echo $V6 $V5 $V7 >> del24.txt
fi
if (( $(echo "0.6 < ${V8} && ${V8} <= 10" | bc -l) )) && [[ "${#V8}" -gt 0 ]]; then
#echo $V6 $V5 $V8
echo $V6 $V5 $V8 >> del48.txt
fi
if (( $(echo "0.6 < ${V9} && ${V9} <= 10" | bc -l) )) && [[ "${#V9}" -gt 0 ]]; then 
#echo $V6 $V5 $V9
echo $V6 $V5 $V9 >> del816.txt
fi
if (( $(echo "0.6 < ${V12} && ${V12} <= 10" | bc -l) )) && [[ "${#V12}" -gt 0 ]]; then
#echo $V6 $V5 $V12
echo $V6 $V5 $V12 >> del1632.txt
fi
fi
fi
done < test.txt
SOD

: << MOD
gmt blockmean del24.txt -R20/29.5/34/39.5 -I0.005=/0.005= > map_24.txt
gmt blockmean del48.txt -R20/29.5/34/39.5 -I0.005=/0.005= > map_48.txt
gmt blockmean del816.txt -R20/29.5/34/39.5 -I0.005=/0.005= > map_816.txt
gmt blockmean del1632.txt -R20/29.5/34/39.5 -I0.005=/0.005= > map_1632.txt

gmt grdmask map_24.txt -R -I0.005=/0.005= -N1/0/0 -Gnew24.nc -S0.25 -V
gmt surface map_24.txt -R -I0.005=/0.005= -GSAegean_grd24.nc -T1.0 -V
gmt grdmask map_48.txt -R -I0.005=/0.005= -N1/0/0 -Gnew48.nc -S0.25 -V
gmt surface map_48.txt -R -I0.005=/0.005= -GSAegean_grd48.nc -T1.0 -V
gmt grdmask map_816.txt -R -I0.005=/0.005= -N1/0/0 -Gnew816.nc -S0.25 -V
gmt surface map_816.txt -R -I0.005=/0.005= -GSAegean_grd816.nc -T1.0 -V
gmt grdmask map_1632.txt -R -I0.005=/0.005= -N1/0/0 -Gnew1632.nc -S0.25 -V
gmt surface map_1632.txt -R -I0.005=/0.005= -GSAegean_grd1632.nc -T1.0 -V
MOD

#: << EOF
PS1=map_24.ps
REGION="20/29.5/34/39.5"
SCALE=1.1i
GRID="SAegean_grd24.nc"
CPT="custom.cpt"

gmt makecpt -Crelief -T-0.24/0.40/0.04 > $CPT

rm -f $PS1

#gmt grdgradient $GRID -A225 -Nt0.5 -G$INT

gmt grdimage $GRID -C$CPT -Inew24.nc -R$REGION -Jm$SCALE -Bf0.5a1/f0.5a1:.:WeSn -V -K -Y1.4 -X2.0 >> $PS1

gmt pscoast -Jm$SCALE -R$REGION -Di -W.3p,0/0/0 -N2p,0/0/0 --DIR_GSHHG=/usr/share/gshhg-gmt-nc4 -V -O -K >> $PS1

gmt psscale -DjRB+w2.5i/0.3i+o2.5/0.5 -R -J -B0.04:P.Delay: -C$CPT -O -V >> $PS1

gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tg
#EOF

PS1=map_48.ps
REGION="20/29.5/34/39.5"
SCALE=1.1i
GRID="SAegean_grd48.nc"
CPT="custom.cpt"

gmt makecpt -Crelief -T-0.24/0.40/0.04 > $CPT

rm -f $PS1

#gmt grdgradient $GRID -A225 -Nt0.5 -G$INT

gmt grdimage $GRID -C$CPT -Inew48.nc -R$REGION -Jm$SCALE -Bf0.5a1/f0.5a1:.:WeSn -V -K -Y1.4 -X2.0 >> $PS1

gmt pscoast -Jm$SCALE -R$REGION -Di -W.3p,0/0/0 -N2p,0/0/0 --DIR_GSHHG=/usr/share/gshhg-gmt-nc4 -V -O -K >> $PS1

gmt psscale -DjRB+w2.5i/0.3i+o2.5/0.5 -R -J -B0.04:P.Delay: -C$CPT -O -V >> $PS1

gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tg

PS1=map_816.ps
REGION="20/29.5/34/39.5"
SCALE=1.1i
GRID="SAegean_grd816.nc"
CPT="custom.cpt"

gmt makecpt -Crelief -T-0.24/0.40/0.04 > $CPT

rm -f $PS1

#gmt grdgradient $GRID -A225 -Nt0.5 -G$INT

gmt grdimage $GRID -C$CPT -Inew48.nc -R$REGION -Jm$SCALE -Bf0.5a1/f0.5a1:.:WeSn -V -K -Y1.4 -X2.0 >> $PS1

gmt pscoast -Jm$SCALE -R$REGION -Di -W.3p,0/0/0 -N2p,0/0/0 --DIR_GSHHG=/usr/share/gshhg-gmt-nc4 -V -O -K >> $PS1

gmt psscale -DjRB+w2.5i/0.3i+o2.5/0.5 -R -J -B0.04:P.Delay: -C$CPT -O -V >> $PS1

gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tg

PS1=map_1632.ps
REGION="20/29.5/34/39.5"
SCALE=1.1i
GRID="SAegean_grd1632.nc"
CPT="custom.cpt"

gmt makecpt -Crelief -T-0.24/0.40/0.04 > $CPT

rm -f $PS1

#gmt grdgradient $GRID -A225 -Nt0.5 -G$INT

gmt grdimage $GRID -C$CPT -Inew1632.nc -R$REGION -Jm$SCALE -Bf0.5a1/f0.5a1:.:WeSn -V -K -Y1.4 -X2.0 >> $PS1

gmt pscoast -Jm$SCALE -R$REGION -Di -W.3p,0/0/0 -N2p,0/0/0 --DIR_GSHHG=/usr/share/gshhg-gmt-nc4 -V -O -K >> $PS1

gmt psscale -DjRB+w2.5i/0.3i+o2.5/0.5 -R -J -B0.04:P.Delay: -C$CPT -O -V >> $PS1

gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tg

