#!/bin/bash

#rm -f index_track.txt nodes_track.txt nodes_track_2-4.txt nodes_track_4-8.txt nodes_track_8-16.txt nodes_track_16-32.txt
#: <<SOD
matlab -nosplash -nodesktop << EOF
fid = fopen('st_delay.txt','r');
A = textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

lonlow = 20.5;
lonhig = 29.5;
latlow = 33.5;
lathig = 39.5;
deplow = 0.0;
dephig = 200.0;
%default num_points_lat = 120, num_points_long = 180, num_points_dep = 6
num_points_lat = 50;
num_points_long = 25;
num_points_dep = 4;
Lat = linspace(latlow,lathig,num_points_lat);
Long = linspace(lonlow,lonhig,num_points_long);
Depth = [135.0000 75.0000 45.0000 15.0000]
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

fid24 = fopen('nodes_24.txt','w');
fid48 = fopen('nodes_48.txt','w');
fid816 = fopen('nodes_816.txt','w');
fid1632 = fopen('nodes_1632.txt','w');


for i = 1:length(A{2})
  %tic
  i
    points = raypath(A{2}(i),A{3}(i),A{5}(i),A{6}(i),A{4}(i));
    index = knnsearch(mdl,points,'K',1);
    bin = bincord(index(:,1),:);
    binf = unique(bin,'rows','stable');
    indexf = unique(index,'rows','stable');
    if (~isnan(A{7}(i)) && A{7}(i) > 0)
      del_ltp = A{7}(i);
      del24 = ones(length(indexf),1)*del_ltp;
      comb24 = [indexf,binf,del24]';
      fprintf(fid24,'%d %f %f %f %f\n',comb24);
    end
    if (~isnan(A{8}(i)) && A{8}(i) > 0)
      del_ltp = A{8}(i);
      del48 = ones(length(indexf),1)*del_ltp;
      comb48 = [indexf,binf,del48]';
      fprintf(fid48,'%d %f %f %f %f\n',comb48);
    end
    if (~isnan(A{9}(i)) && A{9}(i) > 0)
      del_ltp = A{9}(i);
      del816 = ones(length(indexf),1)*del_ltp;
      comb816 = [indexf,binf,del816]';
      fprintf(fid816,'%d %f %f %f %f\n',comb816);
    end
    if (~isnan(A{12}(i)) && A{12}(i) > 0)
      del_ltp = A{12}(i);
      del1632 = ones(length(indexf),1)*del_ltp;
      comb1632 = [indexf,binf,del1632]';
      fprintf(fid1632,'%d %f %f %f %f\n',comb1632);
    end
  %toc
end
fclose(fid24);
fclose(fid48);
fclose(fid816);
fclose(fid1632);
EOF

#SOD

matlab -nosplash -nodesktop << EOF

fid = fopen('nodes_24.txt','r');
A = textscan(fid,'%f %f %f %f %f');
fclose(fid);
comb = [A{1},A{2},A{3},A{4},A{5}];
B = unique(comb(:,1:4),'rows','stable');
C = countmember(B(:,1),comb(:,1));
D = [B,C];
for i=1:length(D)
  a = find(comb(:,1)==D(i,1));
  delay(i) = sum(comb(a,5))/D(i,5);
  i
end
delay =delay';
E = [B,C,delay]';
fid = fopen('nodes_24d.txt','w');
fprintf(fid,'%d %f %f %f %d %f\n',E);
fclose(fid);
EOF

matlab -nosplash -nodesktop << EOF

fid = fopen('nodes_48.txt','r');
A = textscan(fid,'%f %f %f %f %f');
fclose(fid);
comb = [A{1},A{2},A{3},A{4},A{5}];
B = unique(comb(:,1:4),'rows','stable');
C = countmember(B(:,1),comb(:,1));
D = [B,C];
for i=1:length(D)
  a = find(comb(:,1)==D(i,1));
  delay(i) = sum(comb(a,5))/D(i,5);
  i
end
delay =delay';
E = [B,C,delay]';
fid = fopen('nodes_48d.txt','w');
fprintf(fid,'%d %f %f %f %d %f\n',E);
fclose(fid);
EOF

matlab -nosplash -nodesktop << EOF

fid = fopen('nodes_816.txt','r');
A = textscan(fid,'%f %f %f %f %f');
fclose(fid);
comb = [A{1},A{2},A{3},A{4},A{5}];
B = unique(comb(:,1:4),'rows','stable');
C = countmember(B(:,1),comb(:,1));
D = [B,C];
for i=1:length(D)
  a = find(comb(:,1)==D(i,1));
  delay(i) = sum(comb(a,5))/D(i,5);
  i
end
delay =delay';
E = [B,C,delay]';
fid = fopen('nodes_816d.txt','w');
fprintf(fid,'%d %f %f %f %d %f\n',E);
fclose(fid);
EOF

matlab -nosplash -nodesktop << EOF

fid = fopen('nodes_1632.txt','r');
A = textscan(fid,'%f %f %f %f %f');
fclose(fid);
comb = [A{1},A{2},A{3},A{4},A{5}];
B = unique(comb(:,1:4),'rows','stable');
C = countmember(B(:,1),comb(:,1));
D = [B,C];
for i=1:length(D)
  a = find(comb(:,1)==D(i,1));
  delay(i) = sum(comb(a,5))/D(i,5);
  i
end
delay =delay';
E = [B,C,delay]';
fid = fopen('nodes_1632d.txt','w');
fprintf(fid,'%d %f %f %f %d %f\n',E);
fclose(fid);
EOF


