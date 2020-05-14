close all;
clear all;
fid = fopen('st_delay.txt','r');
A = textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

I = find(A{5}==37.07274 | A{5}==37.07);

evlat = A{2}(I);
evlon = A{3}(I);
evdep = A{4}(I);
stlat = A{5}(I);
stlon = A{6}(I);
del24 = A{7}(I);
del48 = A{8}(I);
del816 = A{9}(I);
del812 = A{10}(I);
del1216 = A{11}(I);
del1632 = A{12}(I);
mw =size(I);
dist = zeros(1,mw(1));
az = zeros(1,mw(1));
for i = 1:mw(1)
hdist = distance('gc',A{2}(i),A{3}(i),A{5}(i),A{6}(i))*111.1949;
vdist = A{4}(i);
dist(i) = sqrt(hdist^2 + vdist^2);
az(i) = azimuth('gc',A{5}(i),A{6}(i),A{2}(i),A{3}(i));
end
size(dist)
F = [evlat,evlon,evdep,stlat,stlon,del24,del48,del816,del812,del1216,del1632,dist',az'];
az1 = 121;
az2 = 124;
[lati1,longi1] = reckon(37.0727,25.5230,50/111.1949,az1);
[lati2,longi2] = reckon(37.0727,25.5230,200/111.1949,az1);
[lati3,longi3] = reckon(37.0727,25.5230,400/111.1949,az1);
[latf3,longf3] = reckon(37.0727,25.5230,50/111.1949,az2);
[latf2,longf2] = reckon(37.0727,25.5230,200/111.1949,az2);
[latf1,longf1] = reckon(37.0727,25.5230,400/111.1949,az2);
xv = [25.5230 longi1 longi2 longi3 longf1 longf2 longf3 25.5230];
yv = [37.0727 lati1 lati2 lati3 latf1 latf2 latf3 37.0727];
I = inpolygon(F(:,2),F(:,1),xv,yv);
rw = find(I==1);
figure(1);
scatter(F(rw,12),F(rw,6)); set(gca,'yscale','log'); xlim([30 300]); ylim([0.1 10]); xlabel('Hypocentral Distance (in km)'); ylabel('Delay Time (in s)');
figure(2);
scatter(F(rw,2),F(rw,1)); hold on; plot(xv,yv); hold off; ylabel('Latitude (in deg)'); xlabel('Longitude (in deg)');
figure(3);
scatter3(F(:,2),az,F(:,1)); ylim([90 350]); yticks([90 110 130 150 170 190 210 230 250 270 290 310 330 350]); 
xlabel('Longitude (in deg)'); ylabel('Azimuth (in deg)');
figure(4);
scatter(F(:,2),F(:,1)); xlabel('Longitude (in deg)'); ylabel('Latitude (in deg)');