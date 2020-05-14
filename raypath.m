function [points,thetar,raylen,locs] = raypath(latsrce,lonsrce,latstat,lonstat,depth,num_points)
az = azimuth('gc',latsrce,lonsrce,latstat,lonstat);
odist = distance('gc',latsrce,lonsrce,latstat,lonstat)*111.1949;
z_layer = [0 5 10 15 20 25 30 35 40 50 70 100 150 200 250];
v_layer = [3.08 3.38 3.45 3.45 3.77 3.93 4.10 4.10 4.56 4.59 4.68 4.91 4.91 4.91 4.91];
[val,inf] = closest_value(z_layer,depth);
res = zeros(1,inf);
start = 0;
last = 90;
diff = 1;
j = 0;


while abs(diff) > 0.00001
    med = (start + last)/2;
    theta = pi/180*med;
    res(inf) = (depth-val)*tan(theta);
    for i=inf-1:-1:1
        temp = v_layer(i)*sin(theta)/v_layer(inf);
        res(i) = (z_layer(i+1) - z_layer(i))*temp/sqrt(1-temp^2);
    end
    mdist = sum(res);
    diff = mdist - odist;
    if mdist > odist
        last = med;
    else
        start = med;
    end
    j=j+1;
end


raylen = zeros(1,inf);
thetar = zeros(1,inf);
locs = zeros(inf,2);
raylen(inf) = (depth-val)/cos(theta);
thetar(inf) = theta;
hdist = raylen(inf)*sin(thetar(inf))/111.1949;
locs(inf,:) = reckon('gc',latsrce,lonsrce,hdist,az);


for i = inf-1:-1:1
    temp = v_layer(i)*sin(theta)/v_layer(inf);
    thetar(i) = asin(temp);
    raylen(i) = (z_layer(i+1) - z_layer(i))/sqrt(1-temp^2);
    hdist = raylen(i)*sin(thetar(i))/111.1949;
    locs(i,:) = reckon('gc',locs(i+1,1),locs(i+1,2),hdist,az);
end
raylen = fliplr(raylen);
thetar = fliplr(thetar);
locs = flipud(locs);

lat1 = latsrce;
lon1 = lonsrce; 
k = 1;
points = zeros(num_points*inf,3);
points(1,:) = [latsrce,lonsrce,depth];
for i = 1:inf
    dist = distance('gc',lat1,lon1,locs(i,1),locs(i,2));
    az = azimuth('gc',lat1,lon1,locs(i,1),locs(i,2));
    arlen = dist/num_points;
    pdep = arlen/tan(thetar(i))*111.1949;
    i_st = k+1;
    i_ed = k+num_points;
    points(i_st:i_ed,1:2) = track('gc',[lat1,locs(i,1)],[lon1,locs(i,2)],referenceSphere,'degrees',num_points-1);
    for j = k+1:(k+num_points)
        points(j,3) = depth - pdep * (j-k);
    end
    lat1 = locs(i,1);
    lon1 = locs(i,2);
    depth = points(k+num_points,3);
    k = j;
end

end


function [v, inf] = closest_value(arr, val)
% Returns value and index of arr that is closest to val. If several entries
% are equally close, return the first. Works fine up to machine error (e.g.
% [v, i] = closest_value([4.8, 5], 4.9) will return [5, 2], since in float
% representation 4.9 is strictly closer to 5 than 4.8).
% ===============
% Parameter list:
% ===============
% arr : increasingly ordered array
% val : scalar in R


len = length(arr);
inf = 1;
sup = len;

% Binary search for index
while sup - inf > 1
    med = floor((sup + inf)/2);
    
    % Replace >= here with > to obtain the last index instead of the first.
    if arr(med) >= val 
        sup = med;
    else
        inf = med;
    end
end

% Replace < here with <= to obtain the last index instead of the first.
if sup - inf == 1 && abs(arr(sup) - val) <= abs(arr(inf) - val)
    sup = inf;
end  

v = arr(inf);
end
