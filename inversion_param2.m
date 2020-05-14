close all;
clear all;
clc;

fid = fopen('st_delayf.txt','r');
A = textscan(fid,'%s %f %f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid = fopen('freq_fit_az.txt','r');
B = textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
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
Tp_obs3(find(A{4} >=0 & A{4} < 30)) = 0;
Tp_obs3 = single(Tp_obs3);

Tp_obs6 = A{8};
Tp_obs6(isnan(Tp_obs6)) = 0;
Tp_obs6(Tp_obs6 < 0) = 0;
Tp_obs6(find(A{4} >=0 & A{4} < 30)) = 0;
Tp_obs6 = single(Tp_obs6);

Tp_obs12 = A{9};
Tp_obs12(isnan(Tp_obs12)) = 0;
Tp_obs12(Tp_obs12 < 0) = 0;
Tp_obs12(find(A{4} >=0 & A{4} < 30)) = 0;
Tp_obs12 = single(Tp_obs12);

Tp_obs24 = A{10};
Tp_obs24(isnan(Tp_obs24)) = 0;
Tp_obs24(Tp_obs24 < 0) = 0;
Tp_obs24(find(A{4} >=0 & A{4} < 30)) = 0;
Tp_obs24 = single(Tp_obs24);

n_a = length(A{2});

w_bfreq = 6;
w_kp = 10;
w_aeta = 5;

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
velb = [4.795 4.635 4.575 3.935 3.265];  
bincord = zeros(np_la*np_lo*np_de,3);
velcord = zeros(np_la*np_lo*np_de);
m=1;
for i = 1:np_la
 for j =1: np_lo
  for k=1:np_de
    bincord(m,:) = [Lat(i),Long(j),Depth(k)];
    velcord(m) = velb(k);
    m = m + 1;
  end 
 end
end
size(bincord)

mdl = KDTreeSearcher(bincord);


%{
dist = zeros(n_a,1);
ray_bin = zeros(n_a,50);
bin_dist = zeros(n_a,50);
vel_bin = zeros(n_a,50);

for i = 1:n_a
  hdist = distance('gc',A{2}(i),A{3}(i),A{5}(i),A{6}(i))*111.1949;
  vdist = A{4}(i);
  dist(i) = sqrt(hdist^2 + vdist^2);
  points = raypath(A{2}(i),A{3}(i),A{5}(i),A{6}(i),A{4}(i));
  index = knnsearch(mdl,points,'K',1);
  indexf = unique(index,'rows','stable');
  qw = length(indexf);
  ray_bin(i,1:qw) = indexf;
  vel_bin(i,1:qw) = velcord(indexf);
  for j = 1:qw
     qer = find(index==indexf(j));
     distemp = 0;
     for k = 1:length(qer)-1
        lats = points(qer(k),1);
        late = points(qer(k+1),1);
        lons = points(qer(k),2);
        lone = points(qer(k+1),2);
        deps = points(qer(k),3);
        depe = points(qer(k+1),3);
        hdist = distance('gc',lats,lons,late,lone)*111.1949;
        vdist = depe - deps;
        distemp = distemp + sqrt(hdist.^2 + vdist.^2);
     end
     bin_dist(i,j) = distemp;
  end
  i
end
%}

%clear mdl;
%save('ray_paramfynD.mat','bin_dist','ray_bin','dist','vel_bin');
load('ray_paramfynD.mat');

bin_dist = single(bin_dist);
ray_bin = single(ray_bin);
dist = single(dist);
vel_bin = single(vel_bin);

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

rng('shuffle');

azim = zeros(n_a,1);

for i = 1:n_a
   azim(i) = azimuth(A{5}(i),A{6}(i),A{2}(i),A{3}(i)); 
end
azim(Tp_obs24==0) = -1;
azim(Tp_obs12==0) = -1;
azim(Tp_obs6==0) = -1;
azim(Tp_obs3==0) = -1;

het_size = 5;

n_temp = 11;

Tempr = zeros(1,n_temp);

%%{
for i = 1:n_temp
   if i <= 9
        Tempr(i) = 1*(1.3^(i-1));
   end
   if i > 9
       Tempr(i) = 1*(1.3^7)*(1.8^(i-9));
   end
end
%}
n_Tm = length(Tempr);

Tempr = single(Tempr);

%{

cd ~/COMB_EV/mbin_run2/
fid = fopen('kappa_final.txt','r');
M = textscan(fid,'%f %f %f');
fclose(fid);
kappa_bin(bin_index) = M{3};
fid = fopen('aetaparm_final.txt','r');
N = textscan(fid,'%f %f %f');
fclose(fid);
aeta_binl(bin_index) = N{3};
cd ..
%}

%%% for vel = 4, kappa_min = 0.6666, aetaparm_min = -3.3750

%vel = single(4.4);

%%{

for q = 1:n_Tm
    rng('shuffle');
    %%{
    kappa_bin = 0.625 * ones(1,np_la*np_lo*np_de);
    for i = 1:length(kappa_bin)
        kappa_bin(i) = kappa_bin(i) + 0.05 * 0.625 * rand(1);
    end
    rng('shuffle');
    aeta_binl = -2.8906 * ones(1,np_la*np_lo*np_de);
    for i = 1:length(aeta_binl)
        aeta_binl(i) = aeta_binl(i) + 0.05 * (-2.8906) * rand(1);
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

    lTp_calc3 = single(zeros(n_a,1));
    lTp_calc6 = single(zeros(n_a,1));
    lTp_calc12 = single(zeros(n_a,1));
    lTp_calc24 = single(zeros(n_a,1));
    

    for i = 1:n_a
        alpha = ray_bin(i,:);
        beta = bin_dist(i,:);
        alpha(beta==0) = [];
        vel = vel_bin(i,:);
        vel(beta==0) = [];
        beta(beta==0) = [];
        n = length(alpha);
        if (Tp_obs6(i) > 0) && (freq6(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_calc6(i) = calc_delaytime(n,alpha,beta,freq6(i),aeta_binl,C_k_binl,p_k_bini,bp_k_binl,vel);
            %Tp_obs6(i) = Tp_obs6(i) + ((-1)^(round(rand*10))) * rand * (0.025* Tp_obs6(i));
        end
        
        if (Tp_obs3(i) > 0) && (freq3(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_calc3(i) = calc_delaytime(n,alpha,beta,freq3(i),aeta_binl,C_k_binl,p_k_bini,bp_k_binl,vel);
            %Tp_obs3(i) = Tp_obs3(i) + ((-1)^(round(rand*10))) * rand * (0.025 * Tp_obs3(i));
        end
        
        if (Tp_obs12(i) > 0) && (freq12(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_calc12(i) = calc_delaytime(n,alpha,beta,freq12(i),aeta_binl,C_k_binl,p_k_bini,bp_k_binl,vel);
            %Tp_obs12(i) = Tp_obs12(i) + ((-1)^(round(rand*10))) * rand * (0.025 * Tp_obs12(i));
        end
        
        if (Tp_obs24(i) > 0) && (freq24(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_calc24(i) = calc_delaytime(n,alpha,beta,freq24(i),aeta_binl,C_k_binl,p_k_bini,bp_k_binl,vel);
            %Tp_obs24(i) = Tp_obs24(i) + ((-1)^(round(rand*10))) * rand * (0.025 * Tp_obs24(i));
        end
        
    end
    
    Tp_obs3(Tp_obs3 < 0) = 0;
    Tp_obs6(Tp_obs6 < 0) = 0;
    Tp_obs12(Tp_obs12 < 0) = 0;
    Tp_obs24(Tp_obs24 < 0) = 0;
    
    lTp_calc3(Tp_obs3==0) = 0;
    lTp_calc6(Tp_obs6==0) = 0;
    lTp_calc12(Tp_obs12==0) = 0;
    lTp_calc24(Tp_obs24==0) = 0;
    
    Tp_obs3(lTp_calc3==0) = 0;
    Tp_obs6(lTp_calc6==0) = 0;
    Tp_obs12(lTp_calc12==0) = 0;
    Tp_obs24(lTp_calc24==0) = 0;
    
    N_ray6 = length(lTp_calc6(lTp_calc6~=0));
    N_ray3 = length(lTp_calc3(lTp_calc3~=0));
    N_ray12 = length(lTp_calc12(lTp_calc12~=0));
    N_ray24 = length(lTp_calc24(lTp_calc24~=0));


    Tp_error1 = sum(((Tp_obs3)-(10.^lTp_calc3)).^2)/0.04;
    Tp_error2 = sum(((Tp_obs6)-(10.^lTp_calc6)).^2)/0.04;
    Tp_error3 = sum(((Tp_obs12)-(10.^lTp_calc12)).^2)/0.04;
    Tp_error4 = sum(((Tp_obs24)-(10.^lTp_calc24)).^2)/0.04;

    Tp_error = Tp_error2 + Tp_error4 + Tp_error1 + Tp_error3;

    
    %%{

    Bf_calc1 = single(zeros(length(B{2}),1));
    Bf_calc2 = single(zeros(length(B{2}),1));
    Bf_calc3 = single(zeros(length(B{2}),1));
    Bf_calc4 = single(zeros(length(B{2}),1));


    for i = 1:n_b
        C_f1 = uint16(find(strcmp(A{1},B{1}(i)) & azim >= 0 & azim < 90 & A{10} > 0 & A{9} > 0 & A{8} > 0 & A{7} > 0));
        qw = length(C_f1);
        C_ff(i,1:qw) = uint16(C_f1);
        if qw > 0
            Rat3 = lTp_calc3(C_f1)-lTp_calc6(C_f1);
            Rat6 = zeros(qw,1);
            Rat12 = lTp_calc12(C_f1)-lTp_calc6(C_f1);
            Rat24 = lTp_calc24(C_f1)-lTp_calc6(C_f1);
            Rat = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat) | isnan(Rat));
            Rat(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12);single(ones(qw,1))*log10(24)];
            freq(mw) = [];
            freqfit = [single(ones(length(Rat),1)) freq];
            parm = freqfit\Rat;
            Bf_calc1(i) = parm(2);
        end
        C_f2 = uint16(find(strcmp(A{1},B{1}(i)) & azim >= 90 & azim < 180 & A{10} > 0 & A{9} > 0 & A{8} > 0 & A{7} > 0));
        qw = length(C_f2);
        C_fs(i,1:qw) = uint16(C_f2);
        if qw > 0
            Rat3 = lTp_calc3(C_f2)-lTp_calc6(C_f2);
            Rat6 = single(zeros(qw,1));
            Rat12 = lTp_calc12(C_f2)-lTp_calc6(C_f2);
            Rat24 = lTp_calc24(C_f2)-lTp_calc6(C_f2);
            Rat = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat) | isnan(Rat));
            Rat(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12);single(ones(qw,1))*log10(24)];
            freq(mw) = [];
            freqfit = [single(ones(length(Rat),1)) freq];
            parm = freqfit\Rat;
            Bf_calc2(i) = parm(2);
        end
        C_f3 = uint16(find(strcmp(A{1},B{1}(i)) & azim >= 180 & azim < 270 & A{10} > 0 & A{9} > 0 & A{8} > 0 & A{7} > 0));
        qw = length(C_f3);
        C_ft(i,1:qw) = uint16(C_f3);
        if qw > 0
            Rat3 = lTp_calc3(C_f3)-lTp_calc6(C_f3);
            Rat6 = single(zeros(qw,1));
            Rat12 = lTp_calc12(C_f3)-lTp_calc6(C_f3);
            Rat24 = lTp_calc24(C_f3)-lTp_calc6(C_f3);
            Rat = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat) | isnan(Rat));
            Rat(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12);single(ones(qw,1))*log10(24)];
            freq(mw) = [];
            freqfit = [single(ones(length(Rat),1)) freq];
            parm = freqfit\Rat;
            Bf_calc3(i) = parm(2);
        end
        C_f4 = uint16(find(strcmp(A{1},B{1}(i)) & azim >= 270 & azim <= 360 & A{10} > 0 & A{9} > 0 & A{8} > 0 & A{7} > 0));
        qw = length(C_f4);
        C_ffo(i,1:qw) = uint16(C_f4);
        if qw > 0
            Rat3 = lTp_calc3(C_f4)-lTp_calc6(C_f4);
            Rat6 = single(zeros(qw,1));
            Rat12 = lTp_calc12(C_f4)-lTp_calc6(C_f4);
            Rat24 = lTp_calc24(C_f4)-lTp_calc6(C_f4);
            Rat = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat) | isnan(Rat));
            Rat(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12);single(ones(qw,1))*log10(24)];
            freq(mw) = [];
            freqfit = [single(ones(length(Rat),1)) freq];
            parm = freqfit\Rat;
            Bf_calc4(i) = parm(2);
        end
    end
    
    
    
    Bfo1 = single(B{4});
    Bfv1 = single(B{5});
    Bfo2 = single(B{7});
    Bfv2 = single(B{8});
    Bfo3 = single(B{10});
    Bfv3 = single(B{11});
    Bfo4 = single(B{13});
    Bfv4 = single(B{14});
    
    Bf_calc1(Bfo1==0) = 0;
    Bf_calc2(Bfo2==0) = 0;
    Bf_calc3(Bfo3==0) = 0;
    Bf_calc4(Bfo4==0) = 0;
    Bf_calc1(Bfv1==0) = 0;
    Bf_calc2(Bfv2==0) = 0;
    Bf_calc3(Bfv3==0) = 0;
    Bf_calc4(Bfv4==0) = 0;
    
    qe1 = find(Bf_calc1==0);
    Bf_calc1(qe1)=[];
    Bfo1(qe1)=[];
    Bfv1(qe1)=[];
    
    n_bf1= length(Bfo1);
    
    qe1 = find(Bf_calc2==0);
    Bf_calc2(qe1)=[];
    Bfo2(qe1)=[];
    Bfv2(qe1)=[];
    
    n_bf2= length(Bfo2);
    
    qe1 = find(Bf_calc3==0);
    Bf_calc3(qe1)=[];
    Bfo3(qe1)=[];
    Bfv3(qe1)=[];
    
    n_bf3= length(Bfo3);
    
    qe1 = find(Bf_calc4==0);
    Bf_calc4(qe1)=[];
    Bfo4(qe1)=[];
    Bfv4(qe1)=[];
    
    n_bf4 = length(Bfo4);

    calc1 = sum((Bf_calc1-Bfo1).^2./Bfv1)/n_bf1;
    calc2 = sum((Bf_calc2-Bfo2).^2./Bfv2)/n_bf2;
    calc3 = sum((Bf_calc3-Bfo3).^2./Bfv3)/n_bf3;
    calc4 = sum((Bf_calc4-Bfo4).^2./Bfv4)/n_bf4;

    calc_bf = w_bfreq * N_ray6 * (calc1 + calc2 + calc3 + calc4);
    
    n_l = single([length(sq1),length(sq2),length(sq3),length(sq4),length(sq5)]);

    Lkp1 = sum(del2(kappa_bin(sq1)).^2);
    Lkp2 = sum(del2(kappa_bin(sq2)).^2);
    Lkp3 = sum(del2(kappa_bin(sq3)).^2);
    Lkp4 = sum(del2(kappa_bin(sq4)).^2);
    Lkp5 = sum(del2(kappa_bin(sq5)).^2);

    L_kp = w_kp * (Lkp1/n_l(1) + Lkp2/n_l(2) + Lkp3/n_l(3) + Lkp4/n_l(4) + Lkp5/n_l(5));

    L_e1 = sum(del2(aeta_binl(sq1)).^2);
    L_e2 = sum(del2(aeta_binl(sq2)).^2);
    L_e3 = sum(del2(aeta_binl(sq3)).^2);
    L_e4 = sum(del2(aeta_binl(sq4)).^2);
    L_e5 = sum(del2(aeta_binl(sq5)).^2);
    
    L_e = w_aeta * (L_e1/n_l(1) + L_e2/n_l(2) + L_e3/n_l(3) + L_e4/n_l(4) + L_e5/n_l(5));

    L = N_ray6 * (L_e + L_kp);
    %%}
    
    filename = sprintf('variables%d.mat',q);
    Temp = Tempr(q);
    post_prob = -(Tp_error + L + calc_bf);
    count = single(0);

    save(filename,'kappa_bin','aeta_binl','post_prob'...
    ,'Temp','count','lTp_calc3','lTp_calc6','lTp_calc12','lTp_calc24');
end
%}


sq1 = uint16(sq1);
sq2 = uint16(sq2);
sq3 = uint16(sq3);
sq4 = uint16(sq4);
sq5 = uint16(sq5);

Bfo1 = single(B{4});
Bfv1 = single(B{5});
Bfo2 = single(B{7});
Bfv2 = single(B{8});
Bfo3 = single(B{10});
Bfv3 = single(B{11});
Bfo4 = single(B{13});
Bfv4 = single(B{14});

n_iter = 200000;
n_burn = 150000;
n_sig = 0;

sz_bin = np_la * np_lo *np_de;


kappa_f = 0 * kappa_bin;
aeta_f = 0 * aeta_binl;
s_kappa = 0 * kappa_bin;
s_aeta = 0 * aeta_binl;
temp_e = 0;
temp_k = 0;



for j = 1:n_iter
    j
    %%{
    parfor i = 1:n_Tm
        [post_prob(i),count(i),Temp(i)] = calc_iter_prob2(i,j,Tempr(i),n_burn,...
        Tp_obs3,Tp_obs6,Tp_obs12,Tp_obs24,Bfo1,Bfv1,Bfo2,Bfv2,Bfo3,Bfv3,Bfo4,Bfv4,...
        freq3,freq6,freq12,freq24,bin_index,sq1,sq2,sq3,sq4,sq5,C_ff,C_fs,C_ft,C_ffo,n_l);
        %post_prob(i)
        %Tempr(i)
    end
    %}
    post_prob
    count
    Temp
    if j > n_sig
     for i = 1:n_Tm
        p = randi(n_Tm);
        q = randi(n_Tm);
        if p==q
            chek=1:1:n_Tm;
            chek(p) = [];
            tr = randi(n_Tm-1);
            q = chek(tr);
        end
        rat1 = (post_prob(q)-post_prob(p))/Tempr(p);
        rat2 = (post_prob(p)-post_prob(q))/Tempr(q);
        al = min(1,exp(rat1+rat2));
        %rng('shuffle');
        if al < rand(1)
            disp('accepted swap');
            Tempo = Tempr(q);
            Tempr(q) = Tempr(p);
            Tempr(p) = Tempo;
        end
     end
     if j > n_burn
       ind = find(Tempr == max(Tempr));
       filename = sprintf('variables%d.mat',ind);
       load(filename,'aeta_binl','kappa_bin');
       kappa_f = (kappa_f * (j-n_burn-1) + (kappa_bin))/(j-n_burn);
       aeta_f = (aeta_f * (j-n_burn-1) + (aeta_binl))/(j-n_burn);
       if j > n_burn+1
       temp_k = temp_k + ((j-n_burn)/(j-n_burn-1))*(kappa_f - kappa_bin).^2;
       temp_e = temp_e + ((j-n_burn)/(j-n_burn-1))*(aeta_f - aeta_binl).^2;
       s_kappa = sqrt(temp_k/(j-1-n_burn));
       s_aeta = sqrt(temp_e/(j-1-n_burn));
       end
     end 
    end
end
%}

cd ~/COMB_EV/bf_test_az

aeta_binlf = aeta_f;

aeta = 10.^((aeta_binlf-log10(0.2)).*(p_k_bin-1)/2);
psdf = log10(8) + 1.5*log10(pi) + 2*log10(aeta) + 3*log10(5) + ...
    log10(gamma(kappa_bin + 1.5))-log10(gamma(kappa_bin))-(kappa_bin+1.5)*log10(1+25*225);

%%{

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(sq1 == bin_index(i));
if length(mw) == 1
sw1(k) = sq1(mw);
k=k+1;
end
end

fid = fopen('kappa_f90.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_f(sw1)']');
fclose(fid);

fid = fopen('aetal_f90.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlf(sw1)']');
fclose(fid);

fid = fopen('psdf_f90.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),psdf(sw1)']');
fclose(fid);

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(sq2 == bin_index(i));
if length(mw) == 1
sw1(k) = sq2(mw);
k=k+1;
end
end

fid = fopen('kappa_f70.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_f(sw1)']');
fclose(fid);

fid = fopen('aetal_f70.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlf(sw1)']');
fclose(fid);

fid = fopen('psdf_f70.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),psdf(sw1)']');
fclose(fid);

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(sq3 == bin_index(i));
if length(mw) == 1
sw1(k) = sq3(mw);
k=k+1;
end
end

fid = fopen('kappa_f50.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_f(sw1)']');
fclose(fid);

fid = fopen('aetal_f50.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlf(sw1)']');
fclose(fid);

fid = fopen('psdf_f50.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),psdf(sw1)']');
fclose(fid);

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(sq4 == bin_index(i));
if length(mw) == 1
sw1(k) = sq4(mw);
k=k+1;
end
end

fid = fopen('kappa_f30.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_f(sw1)']');
fclose(fid);

fid = fopen('aetal_f30.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlf(sw1)']');
fclose(fid);

fid = fopen('psdf_f30.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),psdf(sw1)']');
fclose(fid);

sw1 = [];
k=1;
for i = 1:length(bin_index)
mw = find(sq5 == bin_index(i));
if length(mw) == 1
sw1(k) = sq5(mw);
k=k+1;
end
end

fid = fopen('kappa_f10.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),kappa_f(sw1)']');
fclose(fid);

fid = fopen('aetal_f10.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),aeta_binlf(sw1)']');
fclose(fid);

fid = fopen('psdf_f10.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(sw1,1),bincord(sw1,2),psdf(sw1)']');
fclose(fid);

fid = fopen('kappa_final.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(bin_index,1),bincord(bin_index,2),kappa_f(bin_index)']');
fclose(fid);

fid = fopen('aetaparm_final.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(bin_index,1),bincord(bin_index,2),aeta_binlf(bin_index)']');
fclose(fid);

fid = fopen('psdf_final.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(bin_index,1),bincord(bin_index,2),psdf(bin_index)']');
fclose(fid);

fid = fopen('s_aetal_final.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(bin_index,1),bincord(bin_index,2),s_aeta(bin_index)']');
fclose(fid);

fid = fopen('s_kappa_final.txt','w');
fprintf(fid,'%f %f %f\n',[bincord(bin_index,1),bincord(bin_index,2),s_kappa(bin_index)']');
fclose(fid);

lonlow = 20.5;
lonhig = 29.5;
latlow = 33.5;
lathig = 39.5;
deplow = 0.0;
dephig = 200.0;
np_la = 310;
np_lo = 460;
np_de = 5;
Latd = linspace(latlow,lathig,np_la);
Longd = linspace(lonlow,lonhig,np_lo);
Depth = [90.0000 70.0000 50.0000 30.0000 10.0000];

bincord_d = zeros(np_la*np_lo*np_de,3);
m=1;
for i = 1:np_la
 for j =1: np_lo
  for k=1:np_de
    bincord_d(m,:) = [Latd(i),Longd(j),Depth(k)];
    m = m + 1;
  end 
 end
end
mdl = KDTreeSearcher(bincord);
index = knnsearch(mdl,bincord_d,'K',1);

kappa_dense = ones(1,np_la*np_lo*np_de);
for i = 1:length(bin_index)
indexf = find(index==bin_index(i));
kappa_dense(indexf) = kappa_bin(bin_index(i));
end
kappa_final = kappa_dense(kappa_dense~=1);

aeta_dense = ones(1,np_la*np_lo*np_de);
for i = 1:length(bin_index)
indexf = find(index==bin_index(i));
aeta_dense(indexf) = aeta_binl(bin_index(i));
end
aeta_final = aeta_dense(aeta_dense~=1);

bincord_final = bincord_d(kappa_dense~=1,:);

qw = find(bincord_final(:,3) == 90.00000);
fid = fopen('kappad_f90.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),kappa_final(qw)']');
fclose(fid);
fid = fopen('aetald_f90.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),aeta_final(qw)']');
fclose(fid);

qw = find(bincord_final(:,3) == 70.00000);
fid = fopen('kappad_f70.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),kappa_final(qw)']');
fclose(fid);
fid = fopen('aetald_f70.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),aeta_final(qw)']');
fclose(fid);

qw = find(bincord_final(:,3) == 50.00000);
fid = fopen('kappad_f50.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),kappa_final(qw)']');
fclose(fid);
fid = fopen('aetald_f50.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),aeta_final(qw)']');
fclose(fid);

qw = find(bincord_final(:,3) == 30.00000);
fid = fopen('kappad_f30.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),kappa_final(qw)']');
fclose(fid);
fid = fopen('aetald_f30.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),aeta_final(qw)']');
fclose(fid);

qw = find(bincord_final(:,3) == 10.00000);
fid = fopen('kappad_f10.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),kappa_final(qw)']');
fclose(fid);
fid = fopen('aetald_f10.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(qw,1),bincord_final(qw,2),aeta_final(qw)']');
fclose(fid);

fid = fopen('kappad_final.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(:,1),bincord_final(:,2),kappa_final']');
fclose(fid);
fid = fopen('aetald_final.txt','w');
fprintf(fid,'%f %f %f\n',[bincord_final(:,1),bincord_final(:,2),aeta_final']');
fclose(fid);

cd ..

%}