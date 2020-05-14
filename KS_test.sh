#!/bin/bash

. ~/.bashrc

rm -f KS_result.txt

echo -e "St. name"'\t'"P/F48v816"'\t'"816v1632"'\t'"P/F48v1632"'\t'"48v816"'\t''\t'"P/F48v816"'\t'"48v816" >> KS_result.txt

for i in *_freq.txt
do
stnam=${i:0:4}
stnam=${stnam%_}
echo $stnam
matlab -nosplash <<EOF
A = load('$i');
j=1;
k=1;
l=1;
for i=1:length(A)
 if (A(i,1)==6)
  B(j)=A(i,2);
  j=j+1;
 end
 if (A(i,1)==12)
  C(k)=A(i,2);
  k=k+1;
 end
 if (A(i,1)==24)
  D(l)=A(i,2);
  l=l+1;
 end
end
[H4v8,P4v8]=kstest2(B,C,'alpha',0.05);
[H8v16,P8v16]=kstest2(C,D,'alpha',0.05);
[H4v16,P4v16]=kstest2(B,D,'alpha',0.05);
KS=[H4v8 P4v8
H8v16 P8v16
H4v16 P4v16];
dlmwrite("KS_${stnam}.txt",KS,'delimiter','\t','precision',7);
EOF
{
read V1 V2
read M1 M2
read N1 N2 
} < KS_${stnam}.txt
echo -e ${stnam}"\t""\t"$V1"\t""\t"$V2"\t"$M1"\t""\t"$M2"\t"$N1"\t""\t"$N2 >> KS_result.txt
done


