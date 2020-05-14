close all;
clear all;
clc;

fid = fopen('st_delayf.txt','r');
A = textscan(fid,'%s %f %f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid = fopen('freq_fit.txt','r');
B = textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

n_b = length(B{2});

fid = fopen('predom_freq.txt','r');
E = textscan(fid,'%f %f %f %f','Headerlines',1);
fclose(fid);

n_e = length(E{1});

freq3 = E{1};
freq3 = single(freq3);
freq6 = E{2};
freq6 = single(freq6);
freq12 = E{3};
freq12 = single(freq12);
freq24 = E{4};
freq24 = single(freq24);

Tp_obs3 = A{7};
Tp_obs3(isnan(Tp_obs3)) = 0;
Tp_obs3(Tp_obs3 < 0) = 0;
Tp_obs3 = single(Tp_obs3);

Tp_obs6 = A{8};
Tp_obs6(isnan(Tp_obs6)) = 0;
Tp_obs6(Tp_obs6 < 0) = 0;
Tp_obs6 = single(Tp_obs6);

Tp_obs12 = A{9};
Tp_obs12(isnan(Tp_obs12)) = 0;
Tp_obs12(Tp_obs12 < 0) = 0;
Tp_obs12 = single(Tp_obs12);

Tp_obs24 = A{10};
Tp_obs24(isnan(Tp_obs24)) = 0;
Tp_obs24(Tp_obs24 < 0) = 0;
Tp_obs24 = single(Tp_obs24);

n_a = length(A{2});

w_bfreq = 1;
w_kp = 1;
w_aeta = 1;

lonlow = 20.5;
lonhig = 29.5;
latlow = 33.5;
lathig = 39.5;
deplow = 0.0;
dephig = 200.0;
np_la = 31;
np_lo = 46;
np_de = 5;
Lat = linspace(latlow,lathig,np_la);
Long = linspace(lonlow,lonhig,np_lo);
Depth = [90.0000 70.0000 50.0000 30.0000 10.0000];
bincord = zeros(np_la*np_lo*np_de,3);
m=1;
for i = 1:np_la
 for j =1: np_lo
  for k=1:np_de
    bincord(m,:) = [Lat(i),Long(j),Depth(k)];
    m = m + 1;
  end 
 end
end
size(bincord)

pq1 = find(bincord(:,3)==90.0000);
pq2 = find(bincord(:,3)==70.0000);
pq3 = find(bincord(:,3)==50.0000);
pq4 = find(bincord(:,3)==30.0000);
pq5 = find(bincord(:,3)==10.0000);

load('ray_paramfynD.mat');

bin_dist = single(bin_dist);
ray_bin = single(ray_bin);
dist = single(dist);

bin_index = ray_bin(:);
bin_index = unique(bin_index);
bin_index(bin_index==0) = [];
bin_index = single(bin_index);

aq1 = find(bincord(bin_index,3)==90.0000);
aq2 = find(bincord(bin_index,3)==70.0000);
aq3 = find(bincord(bin_index,3)==50.0000);
aq4 = find(bincord(bin_index,3)==30.0000);
aq5 = find(bincord(bin_index,3)==10.0000);

sq1 = bin_index(aq1);
sq2 = bin_index(aq2);
sq3 = bin_index(aq3);
sq4 = bin_index(aq4);
sq5 = bin_index(aq5);

clear aq1 aq2 aq3 aq4 aq5;


kappai = single([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]);
C_ki = single([0.56 1.06 1.56 2.00 2.28 2.31 2.14 1.90 1.68]);
p_ki = single([1.19 1.38 1.56 1.71 1.83 1.91 1.95 1.98 1.99]);
bp_ki = single([0.10 0.17 0.23 0.28 0.31 0.34 0.36 0.37 0.37]);

kappa_bin = 0.8 * ones(1,np_la*np_lo*np_de);
aeta_binl = -6.0 * ones(1,np_la*np_lo*np_de);

kappa_bini = kappa_bin;
aeta_binli = aeta_binl;

%{
xgv = linspace(20.5,29.5,10);
ygv = linspace(33.5,39.5,7);

for j = 1:length(ygv)-1
for i = 1:length(xgv)-1
polyg = [xgv(i),ygv(j);xgv(i+1),ygv(j);xgv(i+1),ygv(j+1);xgv(i),ygv(j+1)];
In1 = inpolygon(bincord(pq1,2),bincord(pq1,1),polyg(:,1),polyg(:,2));
In2 = inpolygon(bincord(pq2,2),bincord(pq2,1),polyg(:,1),polyg(:,2));
In3 = inpolygon(bincord(pq3,2),bincord(pq3,1),polyg(:,1),polyg(:,2));
In4 = inpolygon(bincord(pq4,2),bincord(pq4,1),polyg(:,1),polyg(:,2));
In5 = inpolygon(bincord(pq5,2),bincord(pq5,1),polyg(:,1),polyg(:,2));
flag = (-1)^(i+j);
if flag > 0 
    kappa_bin(pq1(In1)) = 0.75;
    aeta_binl(pq1(In1)) = -1.75;
else
    kappa_bin(pq1(In1)) = 0.25;
    aeta_binl(pq1(In1)) = -5.25;
end
flag = 0-flag;
if flag > 0 
    kappa_bin(pq2(In2)) = 0.75;
    aeta_binl(pq2(In2)) = -1.75;
else
    kappa_bin(pq2(In2)) = 0.25;
    aeta_binl(pq2(In2)) = -5.25;
end
flag = 0-flag;
if flag > 0 
    kappa_bin(pq3(In3)) = 0.75;
    aeta_binl(pq3(In3)) = -1.75;
else
    kappa_bin(pq3(In3)) = 0.25;
    aeta_binl(pq3(In3)) = -5.25;
end
flag = 0-flag;
if flag > 0 
    kappa_bin(pq4(In4)) = 0.75;
    aeta_binl(pq4(In4)) = -1.75;
else
    kappa_bin(pq4(In4)) = 0.25;
    aeta_binl(pq4(In4)) = -5.25;
end
flag = 0-flag;
if flag > 0 
    kappa_bin(pq5(In5)) = 0.75;
    aeta_binl(pq5(In5)) = -1.75;
else
    kappa_bin(pq5(In5)) = 0.25;
    aeta_binl(pq5(In5)) = -5.25;
end
end
end
%}

cd ~/COMB_EV/mbin_run4/
fid = fopen('kappa_final.txt','r');
F = textscan(fid,'%f %f %f');
fclose(fid);
kappa_bin(bin_index) = F{3};
kappa_bini(bin_index) = F{3};
fid = fopen('aetaparm_final.txt','r');
G = textscan(fid,'%f %f %f');
fclose(fid);
aeta_binl(bin_index) = G{3};
aeta_binli(bin_index) = G{3};
cd ..

for i=1:length(pq1)
    kappa_bin(pq1(i)) = kappa_bin(pq1(i)) + 0.10 * kappa_bin(pq1(i))...
        .* sin(2*pi/300 * bincord(pq1(i),2)*111.1949).* ...
        sin(2*pi/300 * bincord(pq1(i),1)*111.1949);
    aeta_binl(pq1(i)) = aeta_binl(pq1(i)) + 0.10 * aeta_binl(pq1(i))...
        .* sin(2*pi/300 * bincord(pq1(i),2)*111.1949).* ...
        sin(2*pi/300 * bincord(pq1(i),1)*111.1949);
end

for i=1:length(pq2)
    kappa_bin(pq2(i)) = kappa_bin(pq2(i)) + 0.10 * kappa_bin(pq2(i))...
        .* cos(2*pi/300 * bincord(pq2(i),2)*111.1949).* ...
        cos(2*pi/300 * bincord(pq2(i),1)*111.1949);
    aeta_binl(pq2(i)) = aeta_binl(pq2(i)) + 0.10 * aeta_binl(pq2(i))...
        .* cos(2*pi/300 * bincord(pq2(i),2)*111.1949).* ...
        cos(2*pi/300 * bincord(pq2(i),1)*111.1949);
end

for i=1:length(pq3)
    kappa_bin(pq3(i)) = kappa_bin(pq3(i)) + 0.10 * kappa_bin(pq3(i))...
        .* sin(2*pi/300 * bincord(pq3(i),2)*111.1949).* ...
        sin(2*pi/300 * bincord(pq3(i),1)*111.1949);
    aeta_binl(pq3(i)) = aeta_binl(pq3(i)) + 0.10 * aeta_binl(pq3(i))...
        .* sin(2*pi/300 * bincord(pq3(i),2)*111.1949).* ...
        sin(2*pi/300 * bincord(pq3(i),1)*111.1949);
end

for i=1:length(pq4)
    kappa_bin(pq4(i)) = kappa_bin(pq4(i)) + 0.10 * kappa_bin(pq4(i))...
        .* cos(2*pi/300 * bincord(pq4(i),2)*111.1949).* ...
        cos(2*pi/300 * bincord(pq4(i),1)*111.1949);
    aeta_binl(pq4(i)) = aeta_binl(pq4(i)) + 0.10 * aeta_binl(pq4(i))...
        .* cos(2*pi/300 * bincord(pq4(i),2)*111.1949).* ...
        cos(2*pi/300 * bincord(pq4(i),1)*111.1949);
end

for i=1:length(pq5)
    kappa_bin(pq5(i)) = kappa_bin(pq5(i)) + 0.10 * kappa_bin(pq5(i))...
        .* sin(2*pi/300 * bincord(pq5(i),2)*111.1949).* ...
        sin(2*pi/300 * bincord(pq5(i),1)*111.1949);
    aeta_binl(pq5(i)) = aeta_binl(pq5(i)) + 0.10 * aeta_binl(pq5(i))...
        .* sin(2*pi/300 * bincord(pq5(i),2)*111.1949).* ...
        sin(2*pi/300 * bincord(pq5(i),1)*111.1949);
end

kappa_binp = kappa_bin - kappa_bini;
aeta_binlp = aeta_binl - aeta_binli;

%%{

cd ~/COMB_EV/synthetic_test_200k/

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(pq1 == bin_index(i));
if length(mw) == 1
sw1(k) = pq1(mw);
k=k+1;
end
end

fid = fopen('kappa_p90.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_binp(sw1)']');
fclose(fid);

fid = fopen('aetal_p90.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlp(sw1)']');
fclose(fid);

fid = fopen('kappa_i90.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_bin(sw1)']');
fclose(fid);

fid = fopen('aetal_i90.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binl(sw1)']');
fclose(fid);

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(pq2 == bin_index(i));
if length(mw) == 1
sw1(k) = pq2(mw);
k=k+1;
end
end

fid = fopen('kappa_p70.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_binp(sw1)']');
fclose(fid);

fid = fopen('aetal_p70.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlp(sw1)']');
fclose(fid);

fid = fopen('kappa_i70.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_bin(sw1)']');
fclose(fid);

fid = fopen('aetal_i70.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binl(sw1)']');
fclose(fid);

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(pq3 == bin_index(i));
if length(mw) == 1
sw1(k) = pq3(mw);
k=k+1;
end
end

fid = fopen('kappa_p50.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_binp(sw1)']');
fclose(fid);

fid = fopen('aetal_p50.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlp(sw1)']');
fclose(fid);

fid = fopen('kappa_i50.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_bin(sw1)']');
fclose(fid);

fid = fopen('aetal_i50.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binl(sw1)']');
fclose(fid);

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(pq4 == bin_index(i));
if length(mw) == 1
sw1(k) = pq4(mw);
k=k+1;
end
end

fid = fopen('kappa_p30.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_binp(sw1)']');
fclose(fid);

fid = fopen('aetal_p30.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlp(sw1)']');
fclose(fid);

fid = fopen('kappa_i30.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_bin(sw1)']');
fclose(fid);

fid = fopen('aetal_i30.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binl(sw1)']');
fclose(fid);

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(pq5 == bin_index(i));
if length(mw) == 1
sw1(k) = pq5(mw);
k=k+1;
end
end

fid = fopen('kappa_p10.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_binp(sw1)']');
fclose(fid);

fid = fopen('aetal_p10.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlp(sw1)']');
fclose(fid);

fid = fopen('kappa_i10.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_bin(sw1)']');
fclose(fid);

fid = fopen('aetal_i10.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binl(sw1)']');
fclose(fid);

fid = fopen('kappa_initial.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(bin_index,1),bincord(bin_index,2),kappa_bin(bin_index)']');
fclose(fid);

fid = fopen('aetal_initial.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(bin_index,1),bincord(bin_index,2),aeta_binl(bin_index)']');
fclose(fid);

cd ..

%}

%%{
    %qw = rand;
    %for i=1:length(F{2})
    %    kappa_bin(F{1}(i)) = 0.1 + 0.8*rat_F(i)*qw;
    %    aeta_binl(F{1}(i)) = -0.6 - 5*rat_F(i)*qw;
    %end
    %{
    kappa_bin = 0.6666 * ones(1,np_la*np_lo*np_de);
    for i = 1:length(kappa_bin)
        kappa_bin(i) = kappa_bin(i) + 0.05 * 0.6666 * rand(1);
    end
    rng('shuffle');
    aeta_binl = -3.3750 * ones(1,np_la*np_lo*np_de);
    for i = 1:length(aeta_binl)
        aeta_binl(i) = aeta_binl(i) + 0.05 * (-3.3750) * rand(1);
    end
    %}
    %filename = sprintf('variables%d.mat',q);
    
    %load(filename,'kappa_bin','aeta_binl','post_prob'...
    %    ,'Temp','count');
    
    kappa_bin = single(kappa_bin);
    aeta_binl = single(aeta_binl);
    
    C_k_bin = interp1(kappai,C_ki,kappa_bin,'spline');
    p_k_bin = interp1(kappai,p_ki,kappa_bin,'spline');
    bp_k_bin = interp1(kappai,bp_ki,kappa_bin,'spline');
    
    C_k_binl = log10(C_k_bin);
    p_k_bini = 1./p_k_bin;
    bp_k_binl = log10(bp_k_bin);

    lTp_syn3 = single(zeros(n_a,1));
    lTp_syn6 = single(zeros(n_a,1));
    lTp_syn12 = single(zeros(n_a,1));
    lTp_syn24 = single(zeros(n_a,1));
    
    Tp_syn3 = single(zeros(n_a,1));
    Tp_syn6 = single(zeros(n_a,1));
    Tp_syn12 = single(zeros(n_a,1));
    Tp_syn24 = single(zeros(n_a,1));
    
    Tp_orig3 = single(zeros(n_a,1));
    Tp_orig6 = single(zeros(n_a,1));
    Tp_orig12 = single(zeros(n_a,1));
    Tp_orig24 = single(zeros(n_a,1));
    

    for i = 1:n_a
        alpha = ray_bin(i,:);
        beta = bin_dist(i,:);
        alpha(beta==0) = [];
        beta(beta==0) = [];
        n = length(alpha);
        if (Tp_obs6(i) > 0) && (freq6(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_syn6(i) = calc_delaytime(n,alpha,beta,freq6(i),aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
            Tp_syn6(i) = 10^lTp_syn6(i);
            Tp_orig6(i) = Tp_syn6(i);
            Tp_syn6(i) = Tp_syn6(i) + ((-1)^i)*0.10*rand*Tp_syn6(i);
        end
        
        if (Tp_obs3(i) > 0) && (freq3(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_syn3(i) = calc_delaytime(n,alpha,beta,freq3(i),aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
            Tp_syn3(i) = 10^lTp_syn3(i);
            Tp_orig3(i) = Tp_syn3(i);
            Tp_syn3(i) = Tp_syn3(i) + ((-1)^i)*0.10*rand*Tp_syn3(i);
        end
        
        if (Tp_obs12(i) > 0) && (freq12(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_syn12(i) = calc_delaytime(n,alpha,beta,freq12(i),aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
            Tp_syn12(i) = 10^lTp_syn12(i);
            Tp_orig12(i) = Tp_syn12(i);
            Tp_syn12(i) = Tp_syn12(i) + ((-1)^i)*0.10*rand*Tp_syn12(i);
        end
        
        if (Tp_obs24(i) > 0) && (freq24(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_syn24(i) = calc_delaytime(n,alpha,beta,freq24(i),aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
            Tp_syn24(i) = 10^lTp_syn24(i);
            Tp_orig24(i) = Tp_syn24(i);
            Tp_syn24(i) = Tp_syn24(i) + ((-1)^i)*0.10*rand*Tp_syn24(i);
        end
        
    end
   
    
    Tp_syn3(lTp_syn3==0) = 0;
    Tp_syn6(lTp_syn6==0) = 0;
    Tp_syn12(lTp_syn12==0) = 0;
    Tp_syn24(lTp_syn24==0) = 0;
   
    
    N_ray6 = length(lTp_syn6(lTp_syn6~=0));
    N_ray3 = length(lTp_syn3(lTp_syn3~=0));
    N_ray12 = length(lTp_syn12(lTp_syn12~=0));
    N_ray24 = length(lTp_syn24(lTp_syn24~=0));

    Bf_syn1 = single(zeros(length(B{2}),1));
    Bf_syn2 = single(zeros(length(B{2}),1));
    Bf_syn3 = single(zeros(length(B{2}),1));
    Bfvd_syn1 = single(zeros(length(B{2}),1));
    Bfvd_syn2 = single(zeros(length(B{2}),1));
    Bfvd_syn3 = single(zeros(length(B{2}),1));

    for i = 1:n_b
        C_f1 = uint16(find(strcmp(A{1},B{1}(i)) & A{4} >= 30 & A{4} < 60));
        qw = length(C_f1);
        C_ff(i,1:qw) = uint16(C_f1);
        if qw > 0
            Rat3 = lTp_syn3(C_f1)-lTp_syn6(C_f1);
            Rat6 = zeros(qw,1);
            Rat12 = lTp_syn12(C_f1)-lTp_syn6(C_f1);
            Rat24 = lTp_syn24(C_f1)-lTp_syn6(C_f1);
            Rat = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat) | isnan(Rat));
            Rat(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);...
                single(ones(qw,1))*log10(12);single(ones(qw,1))*log10(24)];
            freq(mw) = [];
            freqfit = [single(ones(length(Rat),1)) freq];
            parm = freqfit\Rat;
            Bf_syn1(i) = parm(2);
            Bfvd_syn1(i) = sum(((freqfit*parm)-Rat).^2)/length(Rat);
        end
        C_f2 = uint16(find(strcmp(A{1},B{1}(i)) & A{4} >= 60 & A{4} < 90));
        qw = length(C_f2);
        C_fs(i,1:qw) = uint16(C_f2);
        if qw > 0
            Rat3 = lTp_syn3(C_f2)-lTp_syn6(C_f2);
            Rat6 = single(zeros(qw,1));
            Rat12 = lTp_syn12(C_f2)-lTp_syn6(C_f2);
            Rat24 = lTp_syn24(C_f2)-lTp_syn6(C_f2);
            Rat = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat) | isnan(Rat));
            Rat(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12);single(ones(qw,1))*log10(24)];
            freq(mw) = [];
            freqfit = [single(ones(length(Rat),1)) freq];
            parm = freqfit\Rat;
            Bf_syn2(i) = parm(2);
            Bfvd_syn2(i) = sum(((freqfit*parm)-Rat).^2)/length(Rat);
        end
        C_f3 = uint16(find(strcmp(A{1},B{1}(i)) & A{4} >= 90 & A{4} < 250));
        qw = length(C_f3);
        C_ft(i,1:qw) = uint16(C_f3);
        if qw > 0
            Rat3 = lTp_syn3(C_f3)-lTp_syn6(C_f3);
            Rat6 = single(zeros(qw,1));
            Rat12 = lTp_syn12(C_f3)-lTp_syn6(C_f3);
            Rat24 = lTp_syn24(C_f3)-lTp_syn6(C_f3);
            Rat = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat) | isnan(Rat));
            Rat(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12);single(ones(qw,1))*log10(24)];
            freq(mw) = [];
            freqfit = [single(ones(length(Rat),1)) freq];
            parm = freqfit\Rat;
            Bf_syn3(i) = parm(2);
            Bfvd_syn3(i) = sum(((freqfit*parm)-Rat).^2)/length(Rat);
        end
    end

    filename = 'syn_data.mat';
    save(filename,'Tp_syn3','Tp_syn6','Tp_syn12','Tp_syn24','Bf_syn1',...
        'Bfvd_syn1','Bf_syn2','Bfvd_syn2','Bf_syn3','Bfvd_syn3');