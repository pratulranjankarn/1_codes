#!/bin/bash
#rm -rf St_freq
#mkdir St_freq
: <<SOD
rm -f delhyp_2-4.txt delhyp_4-8.txt delhyp_8-16.txt delhyp_16-32.txt delhyp_8-12.txt delhyp_12-16.txt del24.txt del48.txt del816.txt del1621.txt del812.txt del1216.txt
delayfile="st_delay.txt"
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

#~/.bashrc

awk 'NR>=2' $delayfile > temp.txt
while read V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12;
do

if (( $(echo "${V4} >= 40" | bc -l) )); then

#echo $V2 $V3 $V5 $V6 $V4

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

if (( $(echo "${dist} > 50" | bc -l) )) && (( $(echo "${dist} <= 250" | bc -l) )); then

echo $dist $V4

if (( $(echo "${V7} > 0.6" | bc -l) )) && (( $(echo "${V7} <= 10" | bc -l) )) && (( $(echo "${V8} <= 10" | bc -l) )) && (( $(echo "${V8} > 0.6" | bc -l) )); then
echo $V6 $V5 $dist ${V7} >> delhyp_2-4.txt
t_freq24=$(echo "${V7}/${V8}" | bc -l)
echo "3.0" "${t_freq24}" >> ./St_freq/${V1}_freq.txt
if (( $(echo "${V4} >= 30" | bc -l) )) && (( $(echo "${V4} < 50" | bc -l) )); then
echo "3.0" "${t_freq24}" >> ./St_freq/${V1}_1freq.txt
elif (( $(echo "${V4} >= 50" | bc -l) )) && (( $(echo "${V4} < 70" | bc -l) )); then
echo "3.0" "${t_freq24}" >> ./St_freq/${V1}_2freq.txt
elif (( $(echo "${V4} >= 70" | bc -l) )) && (( $(echo "${V4} < 90" | bc -l) )); then
echo "3.0" "${t_freq24}" >> ./St_freq/${V1}_3freq.txt
fi
fi

if (( $(echo "${V8} > 0.6" | bc -l) )) && (( $(echo "${V8} <= 10" | bc -l) )); then
echo $V6 $V5 $dist ${V8} >> delhyp_4-8.txt
#t_freq48=$(echo "${V8}/${V7}" | bc -l)
echo "6.0" "1.0" >> ./St_freq/${V1}_freq.txt
echo "6.0" "1.0" >> ./St_freq/${V1}_1freq.txt
echo "6.0" "1.0" >> ./St_freq/${V1}_2freq.txt
echo "6.0" "1.0" >> ./St_freq/${V1}_3freq.txt
fi

sw=${#V9}
if [[ "${sw}" -gt 0 ]] && (( $(echo "${V9} > 0.6" | bc -l) )) && (( $(echo "${V9} <= 10" | bc -l) )) && (( $(echo "${V8} <= 10" | bc -l) )) && (( $(echo "${V8} > 0.6" | bc -l) )); then
echo $V6 $V5 $dist ${V9} >> delhyp_8-16.txt
t_freq816=$(echo "${V9}/${V8}" | bc -l)
echo "12.0" ${t_freq816} >> ./St_freq/${V1}_freq.txt
if (( $(echo "${V4} >= 40" | bc -l) )) && (( $(echo "${V4} < 60" | bc -l) )); then
echo "12.0" ${t_freq816} >> ./St_freq/${V1}_1freq.txt
elif (( $(echo "${V4} >= 60" | bc -l) )) && (( $(echo "${V4} < 80" | bc -l) )); then
echo "12.0" ${t_freq816} >> ./St_freq/${V1}_2freq.txt
else
echo "12.0" ${t_freq816} >> ./St_freq/${V1}_3freq.txt
fi
fi

sw=${#V12}
if [[ "${sw}" -gt 0 ]] && (( $(echo "${V12} > 0.6" | bc -l) )) && (( $(echo "${V12} <= 10" | bc -l) )) && (( $(echo "${V8} <= 10" | bc -l) )) && (( $(echo "${V8} > 0.6" | bc -l) )); then
echo $V6 $V5 $dist ${V12} >> delhyp_16-32.txt
t_freq1632=$(echo "${V12}/${V8}" | bc -l)
echo "24.0" ${t_freq1632} >> ./St_freq/${V1}_freq.txt
if (( $(echo "${V4} >= 40" | bc -l) )) && (( $(echo "${V4} < 60" | bc -l) )); then
echo "24.0" ${t_freq1632} >> ./St_freq/${V1}_1freq.txt
elif (( $(echo "${V4} >= 60" | bc -l) )) && (( $(echo "${V4} < 80" | bc -l) )); then
echo "24.0" ${t_freq1632} >> ./St_freq/${V1}_2freq.txt
else
echo "24.0" ${t_freq1632} >> ./St_freq/${V1}_3freq.txt
fi
fi

fi

fi

done < temp.txt
SOD

: << ROD

awk '{ print $3 }' delhyp_2-4.txt > test.txt
awk '{ print $4 }' delhyp_2-4.txt > test2.txt
awk '{ print $1 }' delhyp_2-4.txt > testlg.txt
awk '{ print $2 }' delhyp_2-4.txt > testlt.txt

octave << EOF
close all;
clear all;
clc;

fid=fopen('test.txt','r');
A=textscan(fid,"%f");
fclose(fid);
fid=fopen('test2.txt','r');
B=textscan(fid,"%f");
fclose(fid);
fid=fopen('testlg.txt','r');
C=textscan(fid,"%f");
fclose(fid);
fid=fopen('testlt.txt','r');
D=textscan(fid,"%f");
fclose(fid);

pkg load optim

n=size(A{1})(1)
F=[10*ones(n,1),A{1}];

[P,e_var,r,p_var,y_var]=LinearRegression (log10(F),log10(B{1}));
alpha = 10.^(log10(F)*P);
F_it=[A{1},alpha];
diff=log10(B{1})-log10(alpha);
del_ltp = [C{1},D{1},diff]; 
beta = 10.^(log10(alpha) + sqrt(y_var));
gamma = 10.^(log10(alpha) - sqrt(y_var));
B2 = [A{1},beta];
B3 = [A{1},gamma];
dlmwrite("parm_24.txt",P);
dlmwrite("Regg_24.txt",F_it,'delimiter','\t','precision',7);
dlmwrite("Regg_24+.txt",B2,'delimiter','\t','precision',7);
dlmwrite("Regg_24-.txt",B3,'delimiter','\t','precision',7);
dlmwrite("del24.txt",del_ltp,'delimiter','\t','precision',7);
dlmwrite("hist_24.txt",diff);
EOF

awk '{ print $3 }' delhyp_4-8.txt > test.txt
awk '{ print $4 }' delhyp_4-8.txt > test2.txt
awk '{ print $1 }' delhyp_4-8.txt > testlg.txt
awk '{ print $2 }' delhyp_4-8.txt > testlt.txt

octave << EOF
close all;
clear all;
clc;

fid=fopen('test.txt','r');
A=textscan(fid,"%f");
fclose(fid);
fid=fopen('test2.txt','r');
B=textscan(fid,"%f");
fclose(fid);
fid=fopen('testlg.txt','r');
C=textscan(fid,"%f");
fclose(fid);
fid=fopen('testlt.txt','r');
D=textscan(fid,"%f");
fclose(fid);

pkg load optim

n=size(A{1})(1)
F=[10*ones(n,1),A{1}];

[P,e_var,r,p_var,y_var]=LinearRegression (log10(F),log10(B{1}));
alpha = 10.^(log10(F)*P);
F_it=[A{1},alpha];
diff=log10(B{1})-log10(alpha);
del_ltp = [C{1},D{1},diff];
beta = 10.^(log10(alpha) + sqrt(y_var));
gamma = 10.^(log10(alpha) - sqrt(y_var));
B2 = [A{1},beta];
B3 = [A{1},gamma];
dlmwrite("parm_48.txt",P);
dlmwrite("Regg_48.txt",F_it,'delimiter','\t','precision',7);
dlmwrite("Regg_48+.txt",B2,'delimiter','\t','precision',7);
dlmwrite("Regg_48-.txt",B3,'delimiter','\t','precision',7);
dlmwrite("del48.txt",del_ltp,'delimiter','\t','precision',7);
dlmwrite("hist_48.txt",diff);
EOF

awk '{ print $3 }' delhyp_8-16.txt > test.txt
awk '{ print $4 }' delhyp_8-16.txt > test2.txt
awk '{ print $1 }' delhyp_8-16.txt > testlg.txt
awk '{ print $2 }' delhyp_8-16.txt > testlt.txt

octave << EOF
close all;
clear all;
clc;

fid=fopen('test.txt','r');
A=textscan(fid,"%f");
fclose(fid);
fid=fopen('test2.txt','r');
B=textscan(fid,"%f");
fclose(fid);
fid=fopen('testlg.txt','r');
C=textscan(fid,"%f");
fclose(fid);
fid=fopen('testlt.txt','r');
D=textscan(fid,"%f");
fclose(fid);


pkg load optim

n=size(A{1})(1)
F=[10*ones(n,1),A{1}];

[P,e_var,r,p_var,y_var]=LinearRegression (log10(F),log10(B{1}));
alpha = 10.^(log10(F)*P);
F_it=[A{1},alpha];
diff=log10(B{1})-log10(alpha);
del_ltp = [C{1},D{1},diff];
beta = 10.^(log10(alpha) + sqrt(y_var));
gamma = 10.^(log10(alpha) - sqrt(y_var));
B2 = [A{1},beta];
B3 = [A{1},gamma];
dlmwrite("parm_816.txt",P);
dlmwrite("Regg_816.txt",F_it,'delimiter','\t','precision',7);
dlmwrite("Regg_816+.txt",B2,'delimiter','\t','precision',7);
dlmwrite("Regg_816-.txt",B3,'delimiter','\t','precision',7);
dlmwrite("del816.txt",del_ltp,'delimiter','\t','precision',7);
dlmwrite("hist_816.txt",diff);
EOF

awk '{ print $3 }' delhyp_16-32.txt > test.txt
awk '{ print $4 }' delhyp_16-32.txt > test2.txt
awk '{ print $1 }' delhyp_16-32.txt > testlg.txt
awk '{ print $2 }' delhyp_16-32.txt > testlt.txt

octave << EOF
close all;
clear all;
clc;

fid=fopen('test.txt','r');
A=textscan(fid,"%f");
fclose(fid);
fid=fopen('test2.txt','r');
B=textscan(fid,"%f");
fclose(fid);
fid=fopen('testlg.txt','r');
C=textscan(fid,"%f");
fclose(fid);
fid=fopen('testlt.txt','r');
D=textscan(fid,"%f");
fclose(fid);


pkg load optim

n=size(A{1})(1)
F=[10*ones(n,1),A{1}];

[P,e_var,r,p_var,y_var]=LinearRegression (log10(F),log10(B{1}));
alpha = 10.^(log10(F)*P);
F_it=[A{1},alpha];
diff=log10(B{1})-log10(alpha);
del_ltp = [C{1},D{1},diff];
beta = 10.^(log10(alpha) + sqrt(y_var));
gamma = 10.^(log10(alpha) - sqrt(y_var));
B2 = [A{1},beta];
B3 = [A{1},gamma];
dlmwrite("parm_1632.txt",P);
dlmwrite("Regg_1632.txt",F_it,'delimiter','\t','precision',7);
dlmwrite("Regg_1632+.txt",B2,'delimiter','\t','precision',7);
dlmwrite("Regg_1632-.txt",B3,'delimiter','\t','precision',7);
dlmwrite("del1632.txt",del_ltp,'delimiter','\t','precision',7);
dlmwrite("hist_1632.txt",diff);
EOF


ROD

#: << KOD
gmt set MAP_GRID_PEN_PRIMARY thinnest,grey,-

gmt set FONT_LABEL = 12p,Helvetica,black

gmt set FONT_ANNOT_PRIMARY = 12p,Helvetica,black

rm -f plot.ps hist.ps

awk '{ print $3,$4 }' delhyp_2-4.txt > test.txt

{
read W1
read W2
} < parm_24.txt

rm -f plot24.ps
gmt psbasemap -R40/300/1e-1/2e1 -JX5il/3.5il -BWnse -Bxa1pf3g3 -Bya1pf3g3+l'Delay time (s)' -V -K -Y11.8 -X2.0 >> plot.ps
gmt psxy test.txt -R -J -Sc0.2 -W0.5p,black -N -O -V -K >> plot.ps
gmt psxy Regg_24.txt -R -J -W2p,blue -O -V -K >> plot.ps
gmt psxy Regg_24+.txt -R -J -W2p,green -O -V -K >> plot.ps
gmt psxy Regg_24-.txt -R -J -W2p,red -O -V -K >> plot.ps
echo "Intercept="${W1:0:5} | gmt pstext -R -J -DJ1/1 -F+cBR+f12p,Helvetica,black -O -V -K >> plot.ps
echo "Slope="${W2:0:5} | gmt pstext -R -J -DJ1.8/2 -F+cBR+f12p,Helvetica,black -O -V -K >> plot.ps
echo "2-4 Hz" | gmt pstext -R -J -DJ1.8/2.0 -F+cBL+f12p,Helvetica,black -O -V -K >> plot.ps

rm -f hist24.ps
gmt psbasemap -R-1/1/0/1200 -JX5i/3.5i -Bxa0.2f0.1 -Bya200f50+l'Counts' -BWesn -V -K -Y11.8 -X2.0 >> hist.ps
gmt pshistogram hist_24.txt -R -J -W0.05i -F -L0.5p -B -Gblue -O -V -K >> hist.ps
echo "2-4 Hz" | gmt pstext -R -J -DJ1/1.0 -F+cTR+f12p,Helvetica,black -O -V -K >> hist.ps

awk '{ print $3,$4 }' delhyp_4-8.txt > test.txt

{
read W1
read W2
} < parm_48.txt

rm -f plot48.ps
gmt psbasemap -R40/300/1e-1/2e1 -JX5il/3.5il -Bxa1pf3g3 -Bya1pf3g3 -Bwsne -V -O -K -Y0.0 -X14.3 >> plot.ps
gmt psxy test.txt -R -J -Sc0.2 -W0.5p,black -N -O -V -K >> plot.ps
gmt psxy Regg_48.txt -R -J -W2p,blue -O -V -K >> plot.ps
gmt psxy Regg_48+.txt -R -J -W2p,green -O -V -K >> plot.ps
gmt psxy Regg_48-.txt -R -J -W2p,red -O -V -K >> plot.ps
echo "Intercept="${W1:0:5} | gmt pstext -R -J -DJ1/1 -F+cBR+f12p,Helvetica,black -O -V -K >> plot.ps
echo "Slope="${W2:0:5} | gmt pstext -R -J -DJ1.8/2 -F+cBR+f12p,Helvetica,black -O -V -K >> plot.ps
echo "4-8 Hz" | gmt pstext -R -J -DJ1.8/2.0 -F+cBL+f12p,Helvetica,black -O -V -K >> plot.ps

rm -f hist48.ps
gmt psbasemap -R-1/1/0/1200 -JX5i/3.5i -Bxa0.2f0.1 -Bya200f50 -Bwsne -V -O -K -Y0.0 -X14.3 >> hist.ps
gmt pshistogram hist_48.txt -R -J -W0.05i -F -L0.5p -B -Gblue -O -V -K >> hist.ps
echo "4-8 Hz" | gmt pstext -R -J -DJ1/1.0 -F+cTR+f12p,Helvetica,black -O -V -K >> hist.ps

awk '{ print $3,$4 }' delhyp_8-16.txt > test.txt

{
read W1
read W2
} < parm_816.txt

rm -f plot816.ps
gmt psbasemap -R40/300/1e-1/2e1 -JX5il/3.5il -Bxa2f3g3+l'Hypocentral Distance (km)' -Bya1pf3g3+l'Delay time (s)' -BWSne -V -O -K -Y-10.2 -X-14.3 >> plot.ps
gmt psxy test.txt -R -J -Sc0.2 -W0.5p,black -N -O -V -K >> plot.ps
gmt psxy Regg_816.txt -R -J -W2p,blue -O -V -K >> plot.ps
gmt psxy Regg_816+.txt -R -J -W2p,green -O -V -K >> plot.ps
gmt psxy Regg_816-.txt -R -J -W2p,red -O -V -K >> plot.ps
echo "Intercept="${W1:0:5} | gmt pstext -R -J -DJ1/1 -F+cBR+f12p,Helvetica,black -O -V -K >> plot.ps
echo "Slope="${W2:0:5} | gmt pstext -R -J -DJ1.8/2 -F+cBR+f12p,Helvetica,black -O -V -K >> plot.ps
echo "8-16 Hz" | gmt pstext -R -J -DJ1.8/2.0 -F+cBL+f12p,Helvetica,black -O -V -K >> plot.ps

rm -f hist816.ps
gmt psbasemap -R-1/1/0/1200 -JX5i/3.5i -Bxa0.2f0.1+l'Delta log(tp) (s)' -Bya200f50+l'Counts' -BWSne -V -O -K -Y-10.2 -X-14.3 >> hist.ps
gmt pshistogram hist_816.txt -R -J -W0.05i -F -L0.5p -B -Gblue -O -V -K >> hist.ps
echo "8-16 Hz" | gmt pstext -R -J -DJ1/1.0 -F+cTR+f12p,Helvetica,black -O -V -K >> hist.ps

awk '{ print $3,$4 }' delhyp_16-32.txt > test.txt

{
read W1
read W2
} < parm_1632.txt

rm -f plot1632.ps
gmt psbasemap -R40/300/1e-1/2e1 -JX5il/3.5il -Bxa2f3g3+l'Hypocentral Distance (km)' -Bya1pf3g3 -BwSne -V -O -K -Y0.0 -X14.3 >> plot.ps
gmt psxy test.txt -R -J -Sc0.2 -W0.5p,black -N -O -V -K >> plot.ps
gmt psxy Regg_1632.txt -R -J -W2p,blue -O -V -K >> plot.ps
gmt psxy Regg_1632+.txt -R -J -W2p,green -O -V -K >> plot.ps
gmt psxy Regg_1632-.txt -R -J -W2p,red -O -V -K >> plot.ps
echo "Intercept="${W1:0:5} | gmt pstext -R -J -DJ1/1 -F+cBR+f12p,Helvetica,black -O -V -K >> plot.ps
echo "Slope="${W2:0:5} | gmt pstext -R -J -DJ1.8/2 -F+cBR+f12p,Helvetica,black -O -V -K >> plot.ps
echo "16-32 Hz" | gmt pstext -R -J -DJ1.8/2.0 -F+cBL+f12p,Helvetica,black -O -V >> plot.ps

gmt psconvert plot.ps -A+m3c/3c/3c/3c -Tg

gmt psconvert plot.ps -A+m3c/3c/3c/3c -Tjef

rm -f hist1632.ps
gmt psbasemap -R-1/1/0/1200 -JX5i/3.5i -Bxa0.2f0.1+l'Delta log(tp) (s)' -Bya200f50 -BwSen -V -O -K -Y0.0 -X14.3 >> hist.ps
gmt pshistogram hist_1632.txt -R -J -W0.05i -F -L0.5p -B -Gblue -O -V -K >> hist.ps
echo "16-32 Hz" | gmt pstext -R -J -DJ1/1.0 -F+cTR+f12p,Helvetica,black -O -V >> hist.ps
gmt psconvert hist.ps -A+m3c/3c/3c/3c -Tg

gmt psconvert hist.ps -A+m3c/3c/3c/3c -Tjef
KOD

: << COD
rm -f parm_freq_compil.txt
echo -e "Station""\t" "Latitude""\t" "Longitude""\t" "40-60""\t""\t""\t" "Var1""\t""\t""\t" "60-80""\t""\t""\t" "Var2""\t""\t""\t" "80-250""\t""\t""\t" "Var3" >> parm_freq_compil.txt

while read V1 Vlat Vlong
do
echo $V1
chek=`find ./St_freq -name "${V1}_freq.txt"`
qw=${#chek}
if [[ "${qw}" -gt 0 ]]; then
awk '{print $1}' ./St_freq/${V1}_freq.txt > test.txt
awk '{print $2}' ./St_freq/${V1}_freq.txt > test2.txt
octave<<EOF
close all;
clear all;
clc;
fid = fopen('test.txt','r');
A = textscan(fid,"%f\n");
fclose(fid);
fid = fopen('test2.txt','r');
B = textscan(fid,"%f\n");
fclose(fid);

n = size(B{1})(1);
sum24=0;
sum816=0;
sum1632=0;
sum812=0;
sum1216=0;
count24=0;
count816=0;
count1632=0;
count812=0;
count1216=0;
for i=1:n
  if (A{1}(i)==3) 
    sum24+=B{1}(i);
    count24+=1;
  elseif (A{1}(i)==12) 
    sum816+=B{1}(i);
    count816+=1;
  elseif (A{1}(i)==24) 
    sum1632+=B{1}(i);
    count1632+=1;
%  elseif (A{1}(i)==10) 
%    sum812+=B{1}(i);
%    count812+=1;
%  elseif (A{1}(i)==14)
%    sum1216+=B{1}(i);
%    count1216+=1;
  else
    avg48=1.0;
  endif
end
avg24=sum24/count24;
avg816=sum816/count816;
avg1632=sum1632/count1632;
%avg812=sum812/count812;
%avg1216=sum1216/count1216;
q=4;
C1=[3.0 6.0 12.0 24.0]';
D1=[avg24 avg48 avg816 avg1632]';
E1=[C1,D1];
F = [10*ones(n,1),A{1}];

pkg load optim

[P,e_var,r,p_var,y_var] = LinearRegression (log10(F),log10(B{1}));
alpha = 10.^(log10(F)*P);
B5 = [A{1},alpha];
beta = 10.^(log10(alpha) + sqrt(y_var));
gamma = 10.^(log10(alpha) - sqrt(y_var));
B6 = [A{1},beta];
B7 = [A{1},gamma];
my_var = sqrt(sum(y_var.^2)/n);
P_f = [P;my_var];
dlmwrite("ptrg2_${V1}.txt",E1,'delimiter','\t','precision',7);
dlmwrite("ptrg+_${V1}.txt",B6,'delimiter','\t','precision',7);
dlmwrite("ptrg-_${V1}.txt",B7,'delimiter','\t','precision',7);
dlmwrite("parm_${V1}.txt",P_f);
dlmwrite("ptrg_${V1}.txt",B5,'delimiter','\t','precision',7);
EOF

{
read W1
read W2
read W3
} < parm_${V1}.txt
mv parm_${V1}.txt ./St_freq/
rm -f ./St_freq/${V1}freq.ps
gmt psbasemap -R1/40/1e-1/1e2 -JX10il/6il -Bxa1pf3+l'Frequency (Hz)' -Bya1pf3+l'Ratio tp[f Hz]/tp[4-8 Hz]' -BWS+t"Frequency Response of ${V1}" -K > ./St_freq/${V1}freq.ps
gmt psxy ./St_freq/${V1}_freq.txt -R -J -Sc0.3 -W0.3p,grey -O -K >> ./St_freq/${V1}freq.ps
gmt psxy ptrg_${V1}.txt -R -J -W2p,blue -O -K >> ./St_freq/${V1}freq.ps
gmt psxy ptrg+_${V1}.txt -R -J -W2p,green -O -K >> ./St_freq/${V1}freq.ps
gmt psxy ptrg-_${V1}.txt -R -J -W2p,red -O -K >> ./St_freq/${V1}freq.ps
gmt psxy ptrg2_${V1}.txt -R -J -Ss0.3 -W0.4p,black -Gred -O -K >> ./St_freq/${V1}freq.ps
echo "Intercept="${W1:0:5} | gmt pstext -R -J -DJ3/2 -F+cTR -O -V -K >> ./St_freq/${V1}freq.ps
echo "Slope="${W2:0:5} | gmt pstext -R -J -DJ3.8/3 -F+cTR -O -V >> ./St_freq/${V1}freq.ps
gmt psconvert ./St_freq/${V1}freq.ps -A+m3c/3c/3c/3c -Tg

rm -f ptrg_${V1}.txt ptrg2_${V1}.txt ptrg+_${V1}.txt ptrg-_${V1}.txt test.txt test2.txt
rm -f ./St_freq/${V1}freq.ps
else
W1="0"
W2="0"
W3="0"
fi

chek=`find ./St_freq -name "${V1}_1freq.txt"`
qw=${#chek}
if [[ "${qw}" -gt 0 ]]; then
awk '{print $1}' ./St_freq/${V1}_1freq.txt > test.txt
awk '{print $2}' ./St_freq/${V1}_1freq.txt > test2.txt

octave<<EOF
close all;
clear all;
clc;
fid = fopen('test.txt','r');
A = textscan(fid,"%f\n");
fclose(fid);
fid = fopen('test2.txt','r');
B = textscan(fid,"%f\n");
fclose(fid);

n = size(B{1})(1);
sum24=0;
sum816=0;
sum1632=0;
sum812=0;
sum1216=0;
count24=0;
count816=0;
count1632=0;
count812=0;
count1216=0;
for i=1:n
  if (A{1}(i)==3) 
    sum24+=B{1}(i);
    count24+=1;
  elseif (A{1}(i)==12) 
    sum816+=B{1}(i);
    count816+=1;
  elseif (A{1}(i)==24) 
    sum1632+=B{1}(i);
    count1632+=1;
%  elseif (A{1}(i)==10) 
%    sum812+=B{1}(i);
%    count812+=1;
%  elseif (A{1}(i)==14)
%    sum1216+=B{1}(i);
%    count1216+=1;
  else
    avg48=1.0;
  endif
end
avg24=sum24/count24;
avg816=sum816/count816;
avg1632=sum1632/count1632;
%avg812=sum812/count812;
%avg1216=sum1216/count1216;
q=4;
C1=[3.0 6.0 12.0 24.0]';
D1=[avg24 avg48 avg816 avg1632]';
E1=[C1,D1];
F = [10*ones(n,1),A{1}];

pkg load optim

[P,e_var,r,p_var,y_var] = LinearRegression (log10(F),log10(B{1}));
alpha = 10.^(log10(F)*P);
B5 = [A{1},alpha];
beta = 10.^(log10(alpha) + sqrt(y_var));
gamma = 10.^(log10(alpha) - sqrt(y_var));
B6 = [A{1},beta];
B7 = [A{1},gamma];
my_var = sqrt(sum(y_var.^2)/n);
P_f = [P;my_var];
dlmwrite("ptrg2_${V1}_1.txt",E1,'delimiter','\t','precision',7);
dlmwrite("ptrg+_${V1}_1.txt",B6,'delimiter','\t','precision',7);
dlmwrite("ptrg-_${V1}_1.txt",B7,'delimiter','\t','precision',7);
dlmwrite("parm_${V1}_1.txt",P_f);
dlmwrite("ptrg_${V1}_1.txt",B5,'delimiter','\t','precision',7);
EOF

{
read W1
read Wfreq
read Wvar
} < parm_${V1}_1.txt
mv parm_${V1}_1.txt ./St_freq/
rm -f ./St_freq/${V1}freq_1.ps
gmt psbasemap -R1/40/1e-1/1e2 -JX10il/6il -Bxa1pf3+l'Frequency (Hz)' -Bya1pf3+l'Ratio tp[f Hz]/tp[4-8 Hz]' -BWS+t"Frequency Response of ${V1}" -K > ./St_freq/${V1}freq_1.ps
gmt psxy ./St_freq/${V1}_1freq.txt -R -J -Sc0.3 -W0.3p,grey -O -K >> ./St_freq/${V1}freq_1.ps
gmt psxy ptrg_${V1}_1.txt -R -J -W2p,blue -O -K >> ./St_freq/${V1}freq_1.ps
gmt psxy ptrg+_${V1}_1.txt -R -J -W2p,green -O -K >> ./St_freq/${V1}freq_1.ps
gmt psxy ptrg-_${V1}_1.txt -R -J -W2p,red -O -K >> ./St_freq/${V1}freq_1.ps
gmt psxy ptrg2_${V1}_1.txt -R -J -Ss0.3 -W0.4p,black -Gred -O -K >> ./St_freq/${V1}freq_1.ps
echo "Intercept="${W1:0:5} | gmt pstext -R -J -DJ3/2 -F+cTR -O -V -K >> ./St_freq/${V1}freq_1.ps
echo "Slope="${Wfreq:0:5} | gmt pstext -R -J -DJ3.8/3 -F+cTR -O -V >> ./St_freq/${V1}freq_1.ps
gmt psconvert ./St_freq/${V1}freq_1.ps -A+m3c/3c/3c/3c -Tg

rm -f ptrg_${V1}_1.txt ptrg2_${V1}_1.txt ptrg+_${V1}_1.txt ptrg-_${V1}_1.txt test.txt test2.txt
rm -f ./St_freq/${V1}freq_1.ps
else
W1="0"
Wfreq="0"
Wvar="0"
fi

chek=`find ./St_freq -name "${V1}_2freq.txt"`
qw=${#chek}
if [[ "${qw}" -gt 0 ]]; then
awk '{print $1}' ./St_freq/${V1}_2freq.txt > test.txt
awk '{print $2}' ./St_freq/${V1}_2freq.txt > test2.txt

octave<<EOF
close all;
clear all;
clc;
fid = fopen('test.txt','r');
A = textscan(fid,"%f\n");
fclose(fid);
fid = fopen('test2.txt','r');
B = textscan(fid,"%f\n");
fclose(fid);

n = size(B{1})(1);
sum24=0;
sum816=0;
sum1632=0;
sum812=0;
sum1216=0;
count24=0;
count816=0;
count1632=0;
count812=0;
count1216=0;
for i=1:n
  if (A{1}(i)==3) 
    sum24+=B{1}(i);
    count24+=1;
  elseif (A{1}(i)==12) 
    sum816+=B{1}(i);
    count816+=1;
  elseif (A{1}(i)==24) 
    sum1632+=B{1}(i);
    count1632+=1;
%  elseif (A{1}(i)==10) 
%    sum812+=B{1}(i);
%    count812+=1;
%  elseif (A{1}(i)==14)
%    sum1216+=B{1}(i);
%    count1216+=1;
  else
    avg48=1.0;
  endif
end
avg24=sum24/count24;
avg816=sum816/count816;
avg1632=sum1632/count1632;
%avg812=sum812/count812;
%avg1216=sum1216/count1216;
q=4;
C1=[3.0 6.0 12.0 24.0]';
D1=[avg24 avg48 avg816 avg1632]';
E1=[C1,D1];
F = [10*ones(n,1),A{1}];

pkg load optim

[P,e_var,r,p_var,y_var] = LinearRegression (log10(F),log10(B{1}));
alpha = 10.^(log10(F)*P);
B5 = [A{1},alpha];
beta = 10.^(log10(alpha) + sqrt(y_var));
gamma = 10.^(log10(alpha) - sqrt(y_var));
B6 = [A{1},beta];
B7 = [A{1},gamma];
my_var = sqrt(sum(y_var.^2)/n);
P_f = [P;my_var];
dlmwrite("ptrg2_${V1}_2.txt",E1,'delimiter','\t','precision',7);
dlmwrite("ptrg+_${V1}_2.txt",B6,'delimiter','\t','precision',7);
dlmwrite("ptrg-_${V1}_2.txt",B7,'delimiter','\t','precision',7);
dlmwrite("parm_${V1}_2.txt",P_f);
dlmwrite("ptrg_${V1}_2.txt",B5,'delimiter','\t','precision',7);
EOF

{
read W1
read W2freq
read W2var
} < parm_${V1}_2.txt
mv parm_${V1}_2.txt ./St_freq/
rm -f ./St_freq/${V1}freq_2.ps
gmt psbasemap -R1/40/1e-1/1e2 -JX10il/6il -Bxa1pf3+l'Frequency (Hz)' -Bya1pf3+l'Ratio tp[f Hz]/tp[4-8 Hz]' -BWS+t"Frequency Response of ${V1}" -K > ./St_freq/${V1}freq_2.ps
gmt psxy ./St_freq/${V1}_2freq.txt -R -J -Sc0.3 -W0.3p,grey -O -K >> ./St_freq/${V1}freq_2.ps
gmt psxy ptrg_${V1}_2.txt -R -J -W2p,blue -O -K >> ./St_freq/${V1}freq_2.ps
gmt psxy ptrg+_${V1}_2.txt -R -J -W2p,green -O -K >> ./St_freq/${V1}freq_2.ps
gmt psxy ptrg-_${V1}_2.txt -R -J -W2p,red -O -K >> ./St_freq/${V1}freq_2.ps
gmt psxy ptrg2_${V1}_2.txt -R -J -Ss0.3 -W0.4p,black -Gred -O -K >> ./St_freq/${V1}freq_2.ps
echo "Intercept="${W1:0:5} | gmt pstext -R -J -DJ3/2 -F+cTR -O -V -K >> ./St_freq/${V1}freq_2.ps
echo "Slope="${W2freq:0:5} | gmt pstext -R -J -DJ3.8/3 -F+cTR -O -V >> ./St_freq/${V1}freq_2.ps
gmt psconvert ./St_freq/${V1}freq_2.ps -A+m3c/3c/3c/3c -Tg

rm -f ptrg_${V1}_2.txt ptrg2_${V1}_2.txt ptrg+_${V1}_2.txt ptrg-_${V1}_2.txt test.txt test2.txt
rm -f ./St_freq/${V1}freq_2.ps
else
W1="0"
W2freq="0"
W2var="0"
fi

chek=`find ./St_freq -name "${V1}_3freq.txt"`
qw=${#chek}
if [[ "${qw}" -gt 0 ]]; then
awk '{print $1}' ./St_freq/${V1}_3freq.txt > test.txt
awk '{print $2}' ./St_freq/${V1}_3freq.txt > test2.txt

octave<<EOF
close all;
clear all;
clc;
fid = fopen('test.txt','r');
A = textscan(fid,"%f\n");
fclose(fid);
fid = fopen('test2.txt','r');
B = textscan(fid,"%f\n");
fclose(fid);

n = size(B{1})(1);
sum24=0;
sum816=0;
sum1632=0;
sum812=0;
sum1216=0;
count24=0;
count816=0;
count1632=0;
count812=0;
count1216=0;
for i=1:n
  if (A{1}(i)==3) 
    sum24+=B{1}(i);
    count24+=1;
  elseif (A{1}(i)==12) 
    sum816+=B{1}(i);
    count816+=1;
  elseif (A{1}(i)==24) 
    sum1632+=B{1}(i);
    count1632+=1;
%  elseif (A{1}(i)==10) 
%    sum812+=B{1}(i);
%    count812+=1;
%  elseif (A{1}(i)==14)
%    sum1216+=B{1}(i);
%    count1216+=1;
  else
    avg48=1.0;
  endif
end
avg24=sum24/count24;
avg816=sum816/count816;
avg1632=sum1632/count1632;
%avg812=sum812/count812;
%avg1216=sum1216/count1216;
q=4;
C1=[3.0 6.0 12.0 24.0]';
D1=[avg24 avg48 avg816 avg1632]';
E1=[C1,D1];
F = [10*ones(n,1),A{1}];

pkg load optim

[P,e_var,r,p_var,y_var] = LinearRegression (log10(F),log10(B{1}));
alpha = 10.^(log10(F)*P);
B5 = [A{1},alpha];
beta = 10.^(log10(alpha) + sqrt(y_var));
gamma = 10.^(log10(alpha) - sqrt(y_var));
B6 = [A{1},beta];
B7 = [A{1},gamma];
my_var = sqrt(sum(y_var.^2)/n);
P_f = [P;my_var];
dlmwrite("ptrg2_${V1}_3.txt",E1,'delimiter','\t','precision',7);
dlmwrite("ptrg+_${V1}_3.txt",B6,'delimiter','\t','precision',7);
dlmwrite("ptrg-_${V1}_3.txt",B7,'delimiter','\t','precision',7);
dlmwrite("parm_${V1}_3.txt",P_f);
dlmwrite("ptrg_${V1}_3.txt",B5,'delimiter','\t','precision',7);
EOF

{
read W1
read W3freq
read W3var
} < parm_${V1}_3.txt
mv parm_${V1}_3.txt ./St_freq/
rm -f ./St_freq/${V1}freq_3.ps
gmt psbasemap -R1/40/1e-1/1e2 -JX10il/6il -Bxa1pf3+l'Frequency (Hz)' -Bya1pf3+l'Ratio tp[f Hz]/tp[4-8 Hz]' -BWS+t"Frequency Response of ${V1}" -K > ./St_freq/${V1}freq_3.ps
gmt psxy ./St_freq/${V1}_3freq.txt -R -J -Sc0.3 -W0.3p,grey -O -K >> ./St_freq/${V1}freq_3.ps
gmt psxy ptrg_${V1}_3.txt -R -J -W2p,blue -O -K >> ./St_freq/${V1}freq_3.ps
gmt psxy ptrg+_${V1}_3.txt -R -J -W2p,green -O -K >> ./St_freq/${V1}freq_3.ps
gmt psxy ptrg-_${V1}_3.txt -R -J -W2p,red -O -K >> ./St_freq/${V1}freq_3.ps
gmt psxy ptrg2_${V1}_3.txt -R -J -Ss0.3 -W0.4p,black -Gred -O -K >> ./St_freq/${V1}freq_3.ps
echo "Intercept="${W1:0:5} | gmt pstext -R -J -DJ3/2 -F+cTR -O -V -K >> ./St_freq/${V1}freq_3.ps
echo "Slope="${W3freq:0:5} | gmt pstext -R -J -DJ3.8/3 -F+cTR -O -V >> ./St_freq/${V1}freq_3.ps
gmt psconvert ./St_freq/${V1}freq_3.ps -A+m3c/3c/3c/3c -Tg

rm -f ptrg_${V1}_3.txt ptrg2_${V1}_3.txt ptrg+_${V1}_3.txt ptrg-_${V1}_3.txt test.txt test2.txt
rm -f ./St_freq/${V1}freq_3.ps
else
W1="0"
W3freq="0"
W3var="0"
fi
#END
echo -e ${V1}"\t" ${Vlat}"\t" ${Vlong}"\t" ${Wfreq}"\t" ${Wvar}"\t" ${W2freq}"\t" ${W2var}"\t" ${W3freq}"\t" ${W3var} >> parm_freq_compil.txt
done < stlist.txt
#COD
