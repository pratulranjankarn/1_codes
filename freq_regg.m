fid = fopen('st_delay15.txt','r');
A = textscan(fid,'%s %f %f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid = fopen('st_listf.txt','r');
B = textscan(fid,'%f %f %s');
fclose(fid);

n_b = length(B{1});

Tp_obs3 = A{7};
Tp_obs3(isnan(Tp_obs3)) = 1;
Tp_obs3(Tp_obs3 <= 0) = 1;
Tp_obs3(find(A{4} < 35)) = 1;
lTp_obs3 = log10(Tp_obs3);
lTp_obs3 = single(lTp_obs3);

Tp_obs6 = A{8};
Tp_obs6(isnan(Tp_obs6)) = 1;
Tp_obs6(Tp_obs6 <= 0) = 1;
Tp_obs6(find(A{4} < 35)) = 1;
lTp_obs6 = log10(Tp_obs6);
lTp_obs6 = single(lTp_obs6);

Tp_obs12 = A{9};
Tp_obs12(isnan(Tp_obs12)) = 1;
Tp_obs12(Tp_obs12 <= 0) = 1;
Tp_obs12(find(A{4} < 35)) = 1;
lTp_obs12 = log10(Tp_obs12);
lTp_obs12 = single(lTp_obs12);

Tp_obs24 = A{10};
Tp_obs24(isnan(Tp_obs24)) = 1;
Tp_obs24(Tp_obs24 <= 0) = 1;
Tp_obs24(find(A{4} < 35)) = 1;
lTp_obs24 = log10(Tp_obs24);
lTp_obs24 = single(lTp_obs24);

n_a = length(A{2});

Bf_obs1 = single(zeros(length(B{2}),1));
Bfv_obs1 = single(zeros(length(B{2}),1));
nw1 = zeros(n_b,1);
Bf_obs2 = single(zeros(length(B{2}),1));
Bfv_obs2 = single(zeros(length(B{2}),1));
nw2 = zeros(n_b,1);
Bf_obs3 = single(zeros(length(B{2}),1));
Bfv_obs3 = single(zeros(length(B{2}),1));
nw3 = zeros(n_b,1);
Bf_obs4 = single(zeros(length(B{2}),1));
Bfv_obs4 = single(zeros(length(B{2}),1));
nw4 = zeros(n_b,1);
%C_ff = zeros(121,85);

%{

azim = zeros(n_a,1);

for i = 1:n_a
   azim(i) = azimuth(A{5}(i),A{6}(i),A{2}(i),A{3}(i)); 
end
azim(Tp_obs24==1) = -1;
azim(Tp_obs12==1) = -1;
azim(Tp_obs6==1) = -1;
azim(Tp_obs3==1) = -1;
%}
for i = 1:n_b
    C_f1 = uint16(find(strcmp(A{1},B{3}(i)) & A{4} >= 35 & A{4} < 60 & A{9} > 0 & A{8} > 0 & A{7} > 0));
    qw = length(C_f1);
    nw1(i) = qw;
    C_ff(i,1:qw) = uint16(C_f1);
    if (qw > 0) && (nw1(i) > 2)
       Rat3 = lTp_obs3(C_f1)-lTp_obs6(C_f1);
       Rat6 = zeros(qw,1);
       Rat12 = lTp_obs12(C_f1)-lTp_obs6(C_f1);
       %Rat24 = lTp_obs24(C_f1)-lTp_obs6(C_f1);
       Rat = [Rat3;Rat6;Rat12];
       mw = find(isinf(Rat) | isnan(Rat));
       Rat(mw) = [];
       freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12)];
       freq(mw) = [];
       freqfit = [single(ones(length(Rat),1)) freq];
       parm = freqfit\Rat;
       Bf_obs1(i) = parm(2);
       Bfv_obs1(i) = (sum((((freqfit*parm)-Rat).^2)/length(Rat)));
    end
end

%%{

for i = 1:n_b
    C_f1 = uint16(find(strcmp(A{1},B{3}(i)) & A{4} >= 60 & A{4} < 90 & A{9} > 0 & A{8} > 0 & A{7} > 0));
    qw = length(C_f1);
    nw2(i) = qw;
    C_fs(i,1:qw) = uint16(C_f1);
    if (qw > 0) && (nw2(i) > 2)
       Rat3 = lTp_obs3(C_f1)-lTp_obs6(C_f1);
       Rat6 = zeros(qw,1);
       Rat12 = lTp_obs12(C_f1)-lTp_obs6(C_f1);
       %Rat24 = lTp_obs24(C_f1)-lTp_obs6(C_f1);
       Rat = [Rat3;Rat6;Rat12];
       mw = find(isinf(Rat) | isnan(Rat));
       Rat(mw) = [];
       freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12)];
       freq(mw) = [];
       freqfit = [single(ones(length(Rat),1)) freq];
       parm = freqfit\Rat;
       Bf_obs2(i) = parm(2);
       Bfv_obs2(i) = (sum((((freqfit*parm)-Rat).^2)/length(Rat)));
    end
end

for i = 1:n_b
    C_f1 = uint16(find(strcmp(A{1},B{3}(i)) & A{4} >= 90 & A{4} < 250 & A{9} > 0 & A{8} > 0 & A{7} > 0));
    qw = length(C_f1);
    nw3(i) = qw;
    C_ft(i,1:qw) = uint16(C_f1);
    if (qw > 0) && (nw3(i) > 2)
       Rat3 = lTp_obs3(C_f1)-lTp_obs6(C_f1);
       Rat6 = zeros(qw,1);
       Rat12 = lTp_obs12(C_f1)-lTp_obs6(C_f1);
       %Rat24 = lTp_obs24(C_f1)-lTp_obs6(C_f1);
       Rat = [Rat3;Rat6;Rat12];
       mw = find(isinf(Rat) | isnan(Rat));
       Rat(mw) = [];
       freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12)];
       freq(mw) = [];
       freqfit = [single(ones(length(Rat),1)) freq];
       parm = freqfit\Rat;
       Bf_obs3(i) = parm(2);
       Bfv_obs3(i) = (sum((((freqfit*parm)-Rat).^2)/length(Rat)));
    end
end

%{

for i = 1:n_b
    C_f1 = uint16(find(strcmp(A{1},B{3}(i)) & A{4} >= 35 & A{4} < 250 & A{10} > 0 & A{9} > 0 & A{8} > 0 & A{7} > 0));
    qw = length(C_f1);
    nw4(i) = qw;
    C_ff(i,1:qw) = uint16(C_f1);
    if (qw > 0) && (nw4(i) > 2)
       Rat3 = lTp_obs3(C_f1)-lTp_obs6(C_f1);
       Rat6 = zeros(qw,1);
       Rat12 = lTp_obs12(C_f1)-lTp_obs6(C_f1);
       Rat24 = lTp_obs24(C_f1)-lTp_obs6(C_f1);
       Rat = [Rat3;Rat6;Rat12;Rat24];
       mw = find(isinf(Rat) | isnan(Rat));
       Rat(mw) = [];
       freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12);single(ones(qw,1))*log10(24)];
       freq(mw) = [];
       freqfit = [single(ones(length(Rat),1)) freq];
       parm = freqfit\Rat;
       Bf_obs4(i) = parm(2);
       Bfv_obs4(i) = (sum((((freqfit*parm)-Rat).^2)/length(Rat)));
    end
end
%}
alpha = string(B{3});
F = [alpha,B{1},B{2},Bf_obs1,Bfv_obs1,nw1,Bf_obs2,Bfv_obs2,nw2,Bf_obs3,Bfv_obs3,nw3];
%}
%{
alpha = string(B{3});
F = [alpha,B{1},B{2},Bf_obs1,Bfv_obs1,nw1];
%}
fid = fopen('freq_fit15.txt','w');
fprintf(fid,'%s\t %s\t %s\t %s %s %s\t %s %s %s\t %s %s %s\n',F');
fclose(fid);

%{
for i = 1:n_b
    C_f1 = uint16(find(strcmp(A{1},B{3}(i)) & A{4} >= 30 & A{4} < 250 & A{10} > 0 & A{9} > 0 & A{8} > 0 & A{7} > 0));
    qw = length(C_f1);
    nw4(i) = qw;
    C_ff(i,1:qw) = uint16(C_f1);
    if (qw > 0) && (nw4(i) > 2)
       Rat3 = lTp_obs3(C_f1)-lTp_obs6(C_f1);
       Rat6 = zeros(qw,1);
       Rat12 = lTp_obs12(C_f1)-lTp_obs6(C_f1);
       Rat24 = lTp_obs24(C_f1)-lTp_obs6(C_f1);
       Rat = [Rat3;Rat6;Rat12;Rat24];
       mw = find(isinf(Rat) | isnan(Rat));
       Rat(mw) = [];
       freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12);single(ones(qw,1))*log10(24)];
       freq(mw) = [];
       freqfit = [single(ones(length(Rat),1)) freq];
       parm = freqfit\Rat;
       Bf_obs4(i) = parm(2);
       Bfv_obs4(i) = (sum((((freqfit*parm)-Rat).^2)/length(Rat)));
    end
end
%}