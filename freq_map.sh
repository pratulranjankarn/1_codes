#!/bin/bash

export HDF5_DISABLE_VERSION_CHECK=2

rm -f ./g_plus.txt g_comb.txt
rm -f ./g_minus.txt gradient.ps gradient.pdf
: <<EOF
for i in ./St_freq/parm*.txt
do
qw=${i:15:4}
sw=${qw%.}
#grep "${sw}" ./St_freq/KS_result.txt > temp.txt
#read V1 V2 V3 V4 V5 V6 V7 < temp.txt
#echo $V2 $V4 $V6
#if [[ "${V2}" -eq 1 ]] || [[ "${V4}" -eq 1 ]] || [[ "${V6}" -eq 1 ]]; then 
echo $sw 
{
read M1 
read M2
} < $i
echo $M1 $M2
#if (( $(echo "${V2} <= 0" | bc -l) )); then
echo $M2 $sw
mw=`grep "${sw} 3*" stlist.txt` 
#echo $mw $V2 >> g_minus.txt
#fi
#if (( $(echo "${V2} > 0" | bc -l) )); then
#echo $V2 $sw
#mw=`grep "${sw} 3*" stlist.txt` 
echo $mw $M2 >> g_comb.txt
echo $mw $M2
#fi
#fi
done
EOF

awk '{ print $2,$3,$1 }' freq_fit_plot.txt > sc_text.txt

awk '{ print $2,$3,$4 }' freq_fit_plot.txt > triangle.txt

gmt gmtset MAP_FRAME_PEN 3
gmt gmtset MAP_FRAME_WIDTH 0.1
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset PS_PAGE_ORIENTATION landscape
gmt gmtset PS_MEDIA=A4
#gmt gmtset FORMAT_GEO_MAP ddd.mmF

PS1=gradient.ps
REGION="20/29.5/34/39.5"
SCALE=1.1i
GRID="/home/user/Downloads/ETOPO1_Bed_g_gmt4.grd/data"
INT="./topo.int"
CPT="custom.cpt"

gmt set MAP_FRAME_PEN = thickest,black

gmt set MAP_FRAME_TYPE = fancy

gmt set FONT_ANNOT_PRIMARY = 14p,Helvetica,black

rm -f test.txt

gmt makecpt -Cgray -T-8000/0/100 > $CPT
gmt makecpt -Cno_green -T-0.48/0.90/0.04 -Z > comb.cpt

rm -f $PS1

#gmt grdgradient $GRID -A225 -Nt0.5 -G$INT

gmt grdimage $GRID -C$CPT -R$REGION -Jm$SCALE -Bf0.5a1/f0.5a1:.:WeSn -V -K -Y1.4 -X2.0 > $PS1

gmt pscoast -Jm$SCALE -R$REGION -W.3p,0/0/0 -N2p,0/0/0 --DIR_GSHHG=/home/user/GSHHG/gshhg-gmt-2.3.7 -Dh -Lg21/34.5+c37+w100+ar+f -V -O -K >> $PS1

echo "Km" | gmt pstext -R -J -DJ3.35/2.75 -F+cBL+f14p,Helvetica,black -O -V -K >> $PS1

gmt psxy triangle.txt -J -R -O -V -K -W.1 -St.6 -Ccomb.cpt >> $PS1

gmt psscale -DjRB+w2.5i/0.3i+o1.75/0.34 -F+gwhite -R -J -B0.20 -By+l"Bfreq" -Ccomb.cpt -O -V -K >> $PS1

gmt pstext sc_text.txt -R -J -F+f6p,Helvetica-Bold,black -O -V >> $PS1

gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tg

gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tjef

evince $PS1

