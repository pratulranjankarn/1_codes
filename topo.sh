#!/bin/bash

rm -f ray_file.txt seis_data2.txt

#awk '{ print $3,$2,$4 }' st_delay.txt > seis_data.txt

#while read V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12;
#do
#echo ">-Z" >> ray_file.txt
#echo $V3 $V2 >> ray_file.txt
#echo $V6 $V5 >> ray_file.txt
#done < st_delay.txt

#<< SOD
for i in ./loc_Int/*2*.hyp
do
evlalo=`head -n +7 $i | tail -n -1 | awk '{print $12,$10}'`
evlalodep=`head -n +7 $i | tail -n -1 | awk '{print $12,$10,$14}'`
echo $evlalo
tail -n +18 $i | head -n -3 | awk '{ print $1 }' > test.txt
while read V1;
do
stlalo=`grep "$V1 " ./stlist.txt | awk '{print $3,$2}'`
echo $stlalo
echo ">-Z" >> ray_file.txt
echo $evlalo >> ray_file.txt
echo $stlalo >> ray_file.txt
echo $evlalodep >> seis_data2.txt
done < test.txt
done
#SOD
echo "ST"

gmt gmtset MAP_FRAME_PEN 3
gmt gmtset MAP_FRAME_WIDTH 0.1
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset PS_PAGE_ORIENTATION landscape
gmt gmtset PS_MEDIA=A4
#gmt gmtset FORMAT_GEO_MAP ddd.mmF

PS1=topo.ps
REGION="20/29.5/34/39.5"
SCALE=1.1i
GRID="/home/user/Downloads/ETOPO1_Bed_g_gmt4.grd/data"
INT="./topo.int"
CPT="custom.cpt"


rm -f test.txt

gmt makecpt -Cgray -T-8000/0/100 -Z > $CPT

gmt set MAP_FRAME_PEN = thickest,black

gmt set MAP_FRAME_TYPE = fancy

gmt set FONT_ANNOT_PRIMARY = 14p,Helvetica,black

gmt set FONT_LABEL = 14p,Helvetica,black

rm -f $PS1

rm -f topo2.png

echo "SD"

#gmt grdgradient $GRID -A225 -Nt0.5 -G$INTz

gmt grdimage $GRID -C$CPT -R$REGION -Jm$SCALE -Bf0.5a1/f0.5a1:.:WeSn -V -K -Y1.4 -X2.0 > $PS1

gmt psxy ray_file.txt -R -J -O -K  -W.5p,6/79/11 >> $PS1
# 6/79/11

gmt pscoast -Jm$SCALE -R$REGION -W.3p,0/0/0 -N2p,0/0/0 --DIR_GSHHG=/home/user/GSHHG/gshhg-gmt-2.3.7 -Dh -Lg28.5/34.5+c37+w100+ar+f -V -O -K >> $PS1
echo "Km" | gmt pstext -R -J -DJ3.35/2.75 -F+cBR+f14p,Helvetica,black -O -V -K >> $PS1

gmt psxy seis_data2.txt -R -J -O -K  -W.1 -Sc.2 -Cseis2.cpt -h15 >> $PS1

awk '{print $1,$2 }' st_listf.txt > stlist.xy
gmt psxy stlist.xy -R -J -O -K -W.5p,black -St.4 -Ggreen >> $PS1

gmt psscale -R -J -DjLB+w2.5i/0.2i+o1.55/0.32 -F+gwhite -B20 -By+l"Depth(km)" -Cseis2.cpt -O >> $PS1

gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tg

gmt psconvert $PS1 -A+m3c/3c/3c/3c -Tjef

evince $PS1
#gs $PS1
