function [points1,points2,points3] = raypath_par(lat1,lon1,lat2,lon2,dep1,i,num_points,theta)
    points1 = zeros(num_points,1);
    points2 = zeros(num_points,1);
    points3 = zeros(num_points,1);
    dist = distance('gc',lat1,lon1,lat2,lon2);
    az = azimuth('gc',lat1,lon1,lat2,lon2);
    arlen = dist/num_points;
    pdep = arlen/tan(theta)*111.1949;
    i_st = (i-1)*num_points + 1;
    i_end = i*num_points + 1;
    for j = i_st:i_end
        [points1(j),points2(j)] = reckon(lat1,lon1,arlen*(j-i_st),az);
        points3(j) = dep1 - pdep * (j-i_st);
    end

end