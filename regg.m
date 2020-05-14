fid = fopen('st_delay15_pertb.txt','r');
A = textscan(fid,'%s %f %f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);
n_a = length(A{4});
dist = zeros(n_a,1);
for i = 1:n_a
hdist = distance('gc',A{2}(i),A{3}(i),A{5}(i),A{6}(i))*111.1949;
vdist = A{4}(i);
dist(i) = sqrt(hdist^2 + vdist^2);
end
x_ini = dist;
x_ini(find(A{4} < 35)) = 0;
x_ini(find(A{7} <= 0)) = 0;
x_ini(find(dist < 50 | dist > 250)) = 0;


y_ini = A{7};
y_ini(find(A{4} < 35)) = 0;
y_ini(find(A{7} <= 0)) = 0;
y_ini(find(dist < 50 | dist > 250)) = 0;
x_ini(y_ini < 0.3) = 0;
y_ini(y_ini < 0.3) = 0;
y_ini(y_ini == 0) = [];
x_ini(x_ini == 0) = [];

x_ini = log10(x_ini);

y_ini = log10(y_ini);

X = [ones(length(x_ini),1) x_ini];

[b,bint] = regress(y_ini,X);
Y_fit = X*b;
Y_pp = X * bint(:,2);
Y_mm = X * bint(:,1);
Y_pm = X * [bint(1,2); bint(2,1)];
Y_mp = X * [bint(1,1); bint(2,2)];
figure(1);
scatter(x_ini,y_ini); hold on; plot(x_ini,Y_fit,'g'); plot(x_ini,Y_pp,'k');...
    plot(x_ini,Y_mm,'r'); hold off;
d_set = [10.^(x_ini) 10.^(y_ini)]';
d_fitset = [10.^(x_ini) 10.^(Y_fit)]';
b_diff = b - bint(:,1);
parmset = [b b_diff];
parmset = round(parmset,3);
d_diff = y_ini - Y_fit;

fid = fopen('Regg_24.txt','w');
fprintf(fid,'%f %f\n',d_fitset);
fclose(fid);

fid = fopen('parm_24.txt','w');
fprintf(fid,'%f %f\n',parmset');
fclose(fid);

fid = fopen('delhyp_2-4.txt','w');
fprintf(fid,'%f %f\n',d_set);
fclose(fid);

fid = fopen('hist_24.txt','w');
fprintf(fid,'%f\n',d_diff);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_ini = dist;
x_ini(find(A{4} < 35)) = 0;
x_ini(find(A{8} <= 0)) = 0;
x_ini(find(dist < 50 | dist > 250)) = 0;


y_ini = A{8};
y_ini(find(A{4} < 35)) = 0;
y_ini(find(A{8} <= 0)) = 0;
y_ini(find(dist < 50 | dist > 250)) = 0;
x_ini(y_ini < 0.3) = 0;
y_ini(y_ini < 0.3) = 0;
y_ini(y_ini == 0) = [];
x_ini(x_ini == 0) = [];

x_ini = log10(x_ini);

y_ini = log10(y_ini);

X = [ones(length(x_ini),1) x_ini];

[b bint] = regress(y_ini,X);
Y_fit = X*b;
Y_pp = X * bint(:,2);
Y_mm = X * bint(:,1);
Y_pm = X * [bint(1,2); bint(2,1)];
Y_mp = X * [bint(1,1); bint(2,2)];
figure(2);
scatter(x_ini,y_ini); hold on; plot(x_ini,Y_fit,'g'); plot(x_ini,Y_pp,'k');...
    plot(x_ini,Y_mm,'r'); hold off;

d_set = [10.^(x_ini) 10.^(y_ini)]';
d_fitset = [10.^(x_ini) 10.^(Y_fit)]';
b_diff = b - bint(:,1);
parmset = [b b_diff];
parmset = round(parmset,3);
d_diff = y_ini - Y_fit;

fid = fopen('Regg_48.txt','w');
fprintf(fid,'%f %f\n',d_fitset);
fclose(fid);

fid = fopen('parm_48.txt','w');
fprintf(fid,'%f %f\n',parmset');
fclose(fid);

fid = fopen('delhyp_4-8.txt','w');
fprintf(fid,'%f %f\n',d_set);
fclose(fid);

fid = fopen('hist_48.txt','w');
fprintf(fid,'%f\n',d_diff);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_ini = dist;
x_ini(find(A{4} < 35)) = 0;
x_ini(find(A{9} <= 0)) = 0;
x_ini(find(dist < 50 | dist > 250)) = 0;

y_ini = A{9};
y_ini(find(A{4} < 35)) = 0;
y_ini(find(A{9} <= 0)) = 0;
y_ini(find(dist < 50 | dist > 250)) = 0;
x_ini(y_ini < 0.3) = 0;
y_ini(y_ini < 0.3) = 0;
x_ini(isnan(y_ini)) = 0;
y_ini(isnan(y_ini)) = 0;
x_ini(x_ini == 0) = [];
y_ini(y_ini == 0) = [];

x_ini = log10(x_ini);

y_ini = log10(y_ini);

X = [ones(length(x_ini),1) x_ini];

[b bint] = regress(y_ini,X);
Y_fit = X*b;
Y_pp = X * bint(:,2);
Y_mm = X * bint(:,1);
Y_pm = X * [bint(1,2); bint(2,1)];
Y_mp = X * [bint(1,1); bint(2,2)];
figure(3);
scatter(x_ini,y_ini); hold on; plot(x_ini,Y_fit,'g'); plot(x_ini,Y_pp,'k');...
    plot(x_ini,Y_mm,'r'); hold off;

d_set = [10.^(x_ini) 10.^(y_ini)]';
d_fitset = [10.^(x_ini) 10.^(Y_fit)]';
b_diff = b - bint(:,1);
parmset = [b b_diff];
parmset = round(parmset,3);
d_diff = y_ini - Y_fit;

fid = fopen('Regg_816.txt','w');
fprintf(fid,'%f %f\n',d_fitset);
fclose(fid);

fid = fopen('parm_816.txt','w');
fprintf(fid,'%f %f\n',parmset');
fclose(fid);

fid = fopen('delhyp_8-16.txt','w');
fprintf(fid,'%f %f\n',d_set);
fclose(fid);

fid = fopen('hist_816.txt','w');
fprintf(fid,'%f\n',d_diff);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_ini = dist;
x_ini(find(A{4} < 35)) = 0;
x_ini(find(A{10} <= 0)) = 0;
x_ini(find(dist < 50 | dist > 250)) = 0;

y_ini = A{10};
y_ini(find(A{4} < 35)) = 0;
y_ini(find(A{10} <= 0)) = 0;
y_ini(find(dist < 50 | dist > 250)) = 0;
x_ini(y_ini < 0.3) = 0;
y_ini(y_ini < 0.3) = 0;
x_ini(isnan(y_ini)) = 0;
y_ini(isnan(y_ini)) = 0;
x_ini(x_ini == 0) = [];
y_ini(y_ini == 0) = [];

x_ini = log10(x_ini);

y_ini = log10(y_ini);

X = [ones(length(x_ini),1) x_ini];

[b bint] = regress(y_ini,X);
Y_fit = X*b;
Y_pp = X * bint(:,2);
Y_mm = X * bint(:,1);
Y_pm = X * [bint(1,2); bint(2,1)];
Y_mp = X * [bint(1,1); bint(2,2)];
figure(4);
scatter(x_ini,y_ini); hold on; plot(x_ini,Y_fit,'g'); plot(x_ini,Y_pp,'k');...
    plot(x_ini,Y_mm,'r'); hold off;

d_set = [10.^(x_ini) 10.^(y_ini)]';
d_fitset = [10.^(x_ini) 10.^(Y_fit)]';
b_diff = b - bint(:,1);
parmset = [b b_diff];
parmset = round(parmset,3);
d_diff = y_ini - Y_fit;

fid = fopen('Regg_1632.txt','w');
fprintf(fid,'%f %f\n',d_fitset);
fclose(fid);

fid = fopen('parm_1632.txt','w');
fprintf(fid,'%f %f\n',parmset');
fclose(fid);

fid = fopen('delhyp_16-32.txt','w');
fprintf(fid,'%f %f\n',d_set);
fclose(fid);

fid = fopen('hist_1632.txt','w');
fprintf(fid,'%f\n',d_diff);
fclose(fid);

%plot(x_ini,Y_pm,'c'); plot(x_ini,Y_mp,'m');

