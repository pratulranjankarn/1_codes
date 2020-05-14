clc;
close all;
fid = fopen('nodes_24e.txt','r');
A = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);
fid = fopen('nodes_24ep.txt','r');
B = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);

k = 1;

for i=1:length(A{1})
idx = find(B{1}==A{1}(i));
if (~isempty(idx))
diff(k) = A{6}(i) - B{6}(idx);
id(k) = A{1}(i);
lat(k) = A{2}(i);
long(k) = A{3}(i);
dep(k) = A{4}(i);
cnt(k) = A{5}(i);
k = k+1;
end
end

diff = diff';
id = id';
lat=lat';
long=long';
dep=dep';
cnt=cnt';

F = [id lat long dep cnt diff];

fid = fopen('nodes_24diff.txt','w');
fprintf(fid,'%f %f %f %f %f %f\n',F');
fclose(fid);
diff = [];
id = [];
lat = [];
long = [];
dep = [];
cnt = [];

fid = fopen('nodes_48e.txt','r');
A = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);
fid = fopen('nodes_48ep.txt','r');
B = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);

k = 1;

for i=1:length(A{1})
idx = find(B{1}==A{1}(i));
if (~isempty(idx))
diff(k) = A{6}(i) - B{6}(idx);
id(k) = A{1}(i);
lat(k) = A{2}(i);
long(k) = A{3}(i);
dep(k) = A{4}(i);
cnt(k) = A{5}(i);
k = k+1;
end
end

diff = diff';
id = id';
lat=lat';
long=long';
dep=dep';
cnt=cnt';

F = [id lat long dep cnt diff];

fid = fopen('nodes_48diff.txt','w');
fprintf(fid,'%f %f %f %f %f %f\n',F');
fclose(fid);
diff = [];
id = [];
lat = [];
long = [];
dep = [];
cnt = [];

fid = fopen('nodes_816e.txt','r');
A = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);
fid = fopen('nodes_816ep.txt','r');
B = textscan(fid,'%f %f %f %f %f %f');
fclose(fid);

k = 1;

for i=1:length(A{1})
idx = find(B{1}==A{1}(i));
if (~isempty(idx))
diff(k) = A{6}(i) - B{6}(idx);
id(k) = A{1}(i);
lat(k) = A{2}(i);
long(k) = A{3}(i);
dep(k) = A{4}(i);
cnt(k) = A{5}(i);
k = k+1;
end
end

diff = diff';
id = id';
lat=lat';
long=long';
dep=dep';
cnt=cnt';


F = [id lat long dep cnt diff];

fid = fopen('nodes_816diff.txt','w');
fprintf(fid,'%f %f %f %f %f %f\n',F');
fclose(fid);
diff = [];
id = [];
lat = [];
long = [];
dep = [];
cnt = [];

