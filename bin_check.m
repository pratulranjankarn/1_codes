function [bin,index,mdl,points] = bin_check(latsrce,lonsrce,latstat,lonstat,varargin)
points = raypath(latsrce,lonsrce,latstat,lonstat,varargin{1});
%z_layer = [0 5 10 15 20 25 30 35 40 50 70 100 150 200 250];
lonlow = 20.5;
lonhig = 29.5;
latlow = 33.5;
lathig = 39.5;
deplow = 0.0;
dephig = 200.0;
%odep = varargin{1};
num_points_lat = 120;
num_points_long = 180;
num_points_dep = 20;
Lat = linspace(latlow,lathig,num_points_lat);
Long = linspace(lonlow,lonhig,num_points_long);
Depth = linspace(deplow,dephig,num_points_dep);
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
size(bincord);
mdl = KDTreeSearcher(bincord);
index = knnsearch(mdl,points,'K',1);
size(bincord(index(:,1),:));
bin = bincord(index(:,1),:);
end

