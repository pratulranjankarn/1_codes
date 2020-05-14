clc;
close all;
clear all;

cd ../SANT


cd Bf_wav2/

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];
dirinfo(1) = [];
dirinfo(1)=[];

m=1;
subdirinfo = cell(length(dirinfo));
for K = 1 : length(dirinfo)
  thisdir = dirinfo(K).name;
  subdirinfo{K} = dir(fullfile(thisdir, '*cut_E_2-4.sacii'));
  delay12info{K} = dir(fullfile(thisdir, 'delay_12*.txt'));
  delay24info{K} = dir(fullfile(thisdir, 'delay_24*.txt'));
  delay14info{K} = dir(fullfile(thisdir, 'delay_14*.txt'));
  for i = 1:length(subdirinfo{K})
  filepath = strcat(subdirinfo{K}(i).folder,'/',subdirinfo{K}(i).name);
  w_bf(:,m) = waveform(filepath,'sac');
  m = m + 1;
  end
end

cd ..

%{
for i=1:length(subdirinfo)
    temp = load(strcat(subdirinfo{i}.folder,'/',subdirinfo{i}.name));
    delay_bf_12(i) = load(strcat(delay12info{i}.folder,'/',delay12info{i}.name));
    delay_bf_24(i) = load(strcat(delay24info{i}.folder,'/',delay24info{i}.name));
    delay_bf_14(i) = load(strcat(delay14info{i}.folder,'/',delay14info{i}.name));
    if length(temp) < 1701
       x_orig = linspace(0,1700,length(temp));
       x_temp = linspace(0,1700,1701);
       temp2 = interp1(x_orig,temp,x_temp,'spline');
       B(:,i) = temp2;
   else
       B(:,i) = temp;
    end
end
%}


cd Du_Af_wav/

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];
dirinfo(1) = [];
dirinfo(1)=[];

subdirinfo2 = cell(length(dirinfo));
m=1;
for K = 1 : length(dirinfo)
  thisdir = dirinfo(K).name;
  subdirinfo2{K} = dir(fullfile(thisdir, '*cut_E_2-4.sacii'));
  delay12info2{K} = dir(fullfile(thisdir, 'delay_12*.txt'));
  delay24info2{K} = dir(fullfile(thisdir, 'delay_24*.txt'));
  delay14info2{K} = dir(fullfile(thisdir, 'delay_14*.txt'));
  for i = 1:length(subdirinfo2{K})
  filepath = strcat(subdirinfo2{K}(i).folder,'/',subdirinfo2{K}(i).name);
  w_af(:,m) = waveform(filepath,'sac');
  m = m + 1;
  end
end

cd ..

w_f = [w_bf w_af];

c = correlation(w_f);

figure(1);
plot(c,'wig');
title('2-4 Hz');

c = xcorr(c);
c = sort(c);

figure(2);
plot(c,'corr');
title('2-4 Hz');

c = adjusttrig(c,'MIN',1);
figure(3);
plot(c,'sha');
title('2-4 Hz');

c = linkage(c);
figure(4);
plot(c,'den');
title('2-4 Hz');
 
c = cluster(c,0.8);
index = find(c,'CLUST',1);
c1 = subset(c,index);
figure(5);
plot(c1,'wig');
title('2-4 Hz');
%{
for i=1:length(subdirinfo2)
   temp = load(strcat(subdirinfo2{i}.folder,'/',subdirinfo2{i}.name));
   delay_af_12(i) = load(strcat(delay12info2{i}.folder,'/',delay12info2{i}.name));
   delay_af_24(i) = load(strcat(delay24info2{i}.folder,'/',delay24info2{i}.name));
   delay_af_14(i) = load(strcat(delay14info2{i}.folder,'/',delay14info2{i}.name));
   A(1:length(temp),i) = temp;
end

cd ..

c_m = find(B(5,:) == 0);
B(:,c_m) = [];
c_m = find(A(5,:) == 0);
A(:,c_m) = [];
n_1 = size(A);
n_2 = size(B);
%C = cell(n_1,n_2);
%L = cell(n_1,n_2);
for i = 1:n_1(2)
    for j = 1:n_2(2)
        [C{i,j},LAGS{i,j}] = xcorr(A(:,i),B(:,j));
        [R{i,j},P{i,j}] = corrcoef(A(:,i),B(:,j));
        corr_pear(i,j) = abs(R{i,j}(1,2));
        corr_prcnt(i,j) = abs(P{i,j}(1,2));
    end
end

[rs,cls] = find(corr_prcnt >= 0.9);

j = 1;
for i = 1:n_2(2)
    qw = find(cls==i);
    if qw > 0
    delay_corr = rs(qw);
    ind = find(corr_prcnt(delay_corr,i)==max(corr_prcnt(delay_corr,i)));
    ind
    set_ind(j,:) = [i delay_corr(ind)];
    set(j,:) = [delay_bf_12(i) delay_af_12(delay_corr(ind))];
    j = j + 1;
    end
end

[h1 p1] = ttest(set(:,1),set(:,2));
[hk1 pk1] = kstest2(set(:,1),set(:,2));

[hs11 ps11] = kstest(set(:,1));
[hs12 ps12] = kstest(set(:,2));

figure(2)
subplot(2,1,1);
h = scatter(1:1:length(set),set(:,1));
ylabel('Peak Delay time');
title('During Unrest 1-2 Hz');
subplot(2,1,2);
scatter(1:1:length(set),set(:,2));
title('After Unrest 1-2 Hz');
ylabel('Peak Delay time');

j = 1;
for i = 1:n_2(2)
    qw = find(cls==i);
    if qw > 0
    delay_corr = rs(qw);
    ind = find(corr_prcnt(delay_corr,i)==max(corr_prcnt(delay_corr,i)));
    ind
    set_ind(j,:) = [i delay_corr(ind)];
    set(j,:) = [delay_bf_24(i) delay_af_24(delay_corr(ind))];
    j = j + 1;
    end
end

[h2 p2] = ttest(set(:,1),set(:,2));
[hk2 pk2] = kstest2(set(:,1),set(:,2));

[hs21 ps21] = kstest(set(:,1));
[hs22 ps22] = kstest(set(:,2));

figure(3)
subplot(2,1,1);
scatter(1:1:length(set),set(:,1));
title('During Unrest 2-4 Hz');
ylabel('Peak Delay time');
subplot(2,1,2);
scatter(1:1:length(set),set(:,2));
title('After Unrest 2-4 Hz');
ylabel('Peak Delay time');

j = 1;
for i = 1:n_2(2)
    qw = find(cls==i);
    if qw > 0
    delay_corr = rs(qw);
    ind = find(corr_prcnt(delay_corr,i)==max(corr_prcnt(delay_corr,i)));
    ind
    set_ind(j,:) = [i delay_corr(ind)];
    set(j,:) = [delay_bf_14(i) delay_af_14(delay_corr(ind))];
    j = j + 1;
    end
end

[h3 p3] = ttest(set(:,1),set(:,2));
[hk3 pk3] = kstest2(set(:,1),set(:,2));

[hs31 ps31] = kstest(set(:,1));
[hs32 ps32] = kstest(set(:,2));

figure(4)
subplot(2,1,1);
scatter(1:1:length(set),set(:,1));
title('During Unrest 1-4 Hz');
ylabel('Peak Delay time');
subplot(2,1,2);
scatter(1:1:length(set),set(:,2));
title('After Unrest 1-4 Hz');
ylabel('Peak Delay time');
%}
   