#!/bin/bash

#rm -f index_track.txt nodes_track.txt nodes_track_2-4.txt nodes_track_4-8.txt nodes_track_8-16.txt nodes_track_16-32.txt
#
#: <<SOD
matlab -nosplash -nodesktop << EOF
fid = fopen('st_delay15_pertb.txt','r');
A = textscan(fid,'%s %f %f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

lonlow = 20.5;
lonhig = 29.5;
latlow = 33.5;
lathig = 39.5;
deplow = 0.0;
dephig = 200.0;
%default num_points_lat = 120, num_points_long = 180, num_points_dep = 6
num_points_lat = 31;
num_points_long = 46;
num_points_dep = 5;
Lat = linspace(latlow,lathig,num_points_lat);
Long = linspace(lonlow,lonhig,num_points_long);
Depth = [90.0000 70.0000 50.0000 30.0000 10.0000]
bincord = zeros(num_points_lat*num_points_long*num_points_dep,3);
m=1;
for i = 1:num_points_lat
 for j =1: num_points_long
  for k=1:num_points_dep
    bincord(m,:) = [Lat(i),Long(j),Depth(k)];
    m = m + 1;
  end 
 end
end
size(bincord)
mdl = KDTreeSearcher(bincord);

Tp_obs3 = A{7};
Tp_obs3(isnan(Tp_obs3)) = 0;
Tp_obs3(Tp_obs3 < 0) = 0;
Tp_obs3(find(A{4} < 35)) = 0;
Tp_obs3 = single(Tp_obs3);

Tp_obs6 = A{8};
Tp_obs6(isnan(Tp_obs6)) = 0;
Tp_obs6(Tp_obs6 < 0) = 0;
Tp_obs6(find(A{4} < 35)) = 0;
Tp_obs6 = single(Tp_obs6);

Tp_obs12 = A{9};
Tp_obs12(isnan(Tp_obs12)) = 0;
Tp_obs12(Tp_obs12 < 0) = 0;
Tp_obs12(find(A{4} < 35)) = 0;
Tp_obs12 = single(Tp_obs12);

Tp_obs24 = A{10};
Tp_obs24(isnan(Tp_obs24)) = 0;
Tp_obs24(Tp_obs24 < 0) = 0;
Tp_obs24(find(A{4} < 35)) = 0;
Tp_obs24 = single(Tp_obs24);

fid24 = fopen('nodes_24.txt','w');
fid48 = fopen('nodes_48.txt','w');
fid816 = fopen('nodes_816.txt','w');
fid1632 = fopen('nodes_1632.txt','w');

fid = fopen('parm_24.txt','r');
parm24 = textscan(fid,'%f %f');
fclose(fid);
trcpt24 = parm24{1}(1);
slpe24 = parm24{1}(2);

fid = fopen('parm_48.txt','r');
parm48 = textscan(fid,'%f %f');
fclose(fid);
trcpt48 = parm48{1}(1);
slpe48 = parm48{1}(2);

fid = fopen('parm_816.txt','r');
parm816 = textscan(fid,'%f %f');
fclose(fid);
trcpt816 = parm816{1}(1);
slpe816 = parm816{1}(2);

fid = fopen('parm_1632.txt','r');
parm1632 = textscan(fid,'%f %f\n');
fclose(fid);
trcpt1632 = parm1632{1}(1);
slpe1632 = parm1632{1}(2);


for i = 1:length(A{2})
  %tic
  i
  hdist = distance('gc',A{2}(i),A{3}(i),A{5}(i),A{6}(i))*111.1949;
  vdist = A{4}(i);
  dist = sqrt(hdist^2 + vdist^2);
  if (dist >= 50 && dist <= 250)
    points = raypath(A{2}(i),A{3}(i),A{5}(i),A{6}(i),A{4}(i),1000);
    [rw,cw] = find(isnan(points));
    points(rw,:) = [];
    index = knnsearch(mdl,points,'K',1);
    bin = bincord(index(:,1),:);
    binf = unique(bin,'rows','stable');
    indexf = unique(index,'rows','stable');
    if (Tp_obs3(i) > 0.3)
      del_ltp = log10(Tp_obs3(i)) - (trcpt24 + slpe24*log10(dist));
      del24 = ones(length(indexf),1)*del_ltp;
      w24 = ones(length(indexf),1)*length(indexf);
      comb24 = [indexf,binf,w24,del24]';
      fprintf(fid24,'%d %f %f %f %d %f\n',comb24);
    end
    if (Tp_obs6(i) > 0.3)
      del_ltp = log10(Tp_obs6(i)) - (trcpt48 + slpe48*log10(dist));
      del48 = ones(length(indexf),1)*del_ltp;
      w48 = ones(length(indexf),1)*length(indexf);
      comb48 = [indexf,binf,w48,del48]';
      fprintf(fid48,'%d %f %f %f %d %f\n',comb48);
    end
    if (Tp_obs12(i) > 0.3)
      del_ltp = log10(Tp_obs12(i)) - (trcpt816 + slpe816*log10(dist));
      del816 = ones(length(indexf),1)*del_ltp;
      w816 = ones(length(indexf),1)*length(indexf);
      comb816 = [indexf,binf,w816,del816]';
      fprintf(fid816,'%d %f %f %f %d %f\n',comb816);
    end
    if (Tp_obs24(i) > 0.3)
      del_ltp = log10(Tp_obs24(i)) - (trcpt1632 + slpe1632*log10(dist));
      del1632 = ones(length(indexf),1)*del_ltp;
      w1632 = ones(length(indexf),1)*length(indexf);
      comb1632 = [indexf,binf,w1632,del1632]';
      fprintf(fid1632,'%d %f %f %f %d %f\n',comb1632);
    end
  end
  %toc
end
fclose(fid24);
fclose(fid48);
fclose(fid816);
fclose(fid1632);
EOF

SOD

matlab -nosplash -nodesktop << EOF

fid = fopen('nodes_24.txt','r');
A = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);
comb = [A{1},A{2},A{3},A{4},A{5},A{6}];
B = unique(comb(:,1:4),'rows','stable');
C = countmember(B(:,1),comb(:,1));
D = [B,C];
for i=1:length(D)
  a = find(comb(:,1)==D(i,1));
  %chek = sum(comb(a,6)./comb(a,5))/C(i);
  chek = min(comb(a,6));
  %chek = min(comb(a,6)./comb(a,5));
  %chek = log10(sum((10.^comb(a,6))./comb(a,5))/C(i));
  if length(a) == C(i)
    disp('equal')
  end
  delay(i) = chek;
  i
end
delay =delay';
E = [B,C,delay]';
fid = fopen('nodes_24ep.txt','w');
fprintf(fid,'%d %f %f %f %d %f\n',E);
fclose(fid);
EOF

matlab -nosplash -nodesktop << EOF

fid = fopen('nodes_48.txt','r');
A = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);
comb = [A{1},A{2},A{3},A{4},A{5},A{6}];
B = unique(comb(:,1:4),'rows','stable');
C = countmember(B(:,1),comb(:,1));
D = [B,C];
for i=1:length(D)
  a = find(comb(:,1)==D(i,1));
  %chek = sum(comb(a,6)./comb(a,5))/length(a);
  chek = min(comb(a,6));
  %chek = min(comb(a,6)./comb(a,5));
  %chek = log10(sum((10.^comb(a,6))./comb(a,5))/C(i));
  if length(a) == C(i)
    disp('equal')
  end
  delay(i) = chek;
  i
end
delay =delay';
E = [B,C,delay]';
fid = fopen('nodes_48ep.txt','w');
fprintf(fid,'%d %f %f %f %d %f\n',E);
fclose(fid);
EOF

matlab -nosplash -nodesktop << EOF

fid = fopen('nodes_816.txt','r');
A = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);
comb = [A{1},A{2},A{3},A{4},A{5},A{6}];
B = unique(comb(:,1:4),'rows','stable');
C = countmember(B(:,1),comb(:,1));
D = [B,C];
for i=1:length(D)
  a = find(comb(:,1)==D(i,1));
  %chek = sum(comb(a,6)./comb(a,5))/length(a);
  chek = min(comb(a,6));
  %chek = min(comb(a,6)./comb(a,5));
  %chek = log10(sum((10.^comb(a,6))./comb(a,5))/C(i));
  if length(a) == C(i)
    disp('equal')
  end
  delay(i) = chek;
  i
end
delay =delay';
E = [B,C,delay]';
fid = fopen('nodes_816ep.txt','w');
fprintf(fid,'%d %f %f %f %d %f\n',E);
fclose(fid);
EOF

matlab -nosplash -nodesktop << EOF

fid = fopen('nodes_1632.txt','r');
A = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);
comb = [A{1},A{2},A{3},A{4},A{5},A{6}];
B = unique(comb(:,1:4),'rows','stable');
C = countmember(B(:,1),comb(:,1));
D = [B,C];
for i=1:length(D)
  a = find(comb(:,1)==D(i,1));
  %chek = sum(comb(a,6)./comb(a,5))/length(a);
  chek = min(comb(a,6));
  %chek = min(comb(a,6)./comb(a,5));
  %chek = log10(sum((10.^comb(a,6))./comb(a,5))/C(i));
  if length(a) == C(i)
    disp('equal')
  end
  delay(i) = chek;
  i
end
delay =delay';
E = [B,C,delay]';
fid = fopen('nodes_1632ep.txt','w');
fprintf(fid,'%d %f %f %f %d %f\n',E);
fclose(fid);
EOF


