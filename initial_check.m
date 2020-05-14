aeta_binli = linspace(-3.8281,-3.7344,9);
kappai = linspace(0.45,0.5,9);

kappi = single([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]);
C_ki = single([0.56 1.06 1.56 2.00 2.28 2.31 2.14 1.90 1.68]);
p_ki = single([1.19 1.38 1.56 1.71 1.83 1.91 1.95 1.98 1.99]);
bp_ki = single([0.10 0.17 0.23 0.28 0.31 0.34 0.36 0.37 0.37]);

k=1;
parm_set = zeros(length(kappai)*length(aeta_binli),2); 

for i = 1:length(kappai)
    for j = 1:length(aeta_binli)
        parm_set(k,:) = [kappai(i),aeta_binli(j)];
        k = k + 1;
    end
end

n_set = length(kappai)*length(aeta_binli);

load('ray_paramcrsD_30.mat');

%fid = fopen('predom_freq.txt','r');
%E = textscan(fid,'%f %f %f %f','Headerlines',1);
%fclose(fid);

%n_e = length(E{1});

freq3 = 3 * ones(length(A{4}),1);
freq3 = single(freq3);
freq6 = 6 * ones(length(A{4}),1);
freq6 = single(freq6);
freq12 = 12 * ones(length(A{4}),1);
freq12 = single(freq12);
freq24 = 24 * ones(length(A{4}),1);
freq24 = single(freq24);

n_a = length(Tp_obs3);


vel_bin = 0 * single(vel_bin) + 4;

for q = 1:n_set
    q
    rng('shuffle');
    %qw = rand;
    %for i=1:length(F{2})
    %    kappa_bin(F{1}(i)) = 0.1 + 0.8*rat_F(i)*qw;
    %    aeta_binl(F{1}(i)) = -0.6 - 5*rat_F(i)*qw;
    %end
    kappa_bin = parm_set(q,1)*ones(1,np_la*np_lo*np_de);
    aeta_binl = parm_set(q,2)*ones(1,np_la*np_lo*np_de);

    %filename = sprintf('variables%d.mat',q);
    
    %load(filename,'kappa_bin','aeta_binl','post_prob'...
    %    ,'Temp','count');
    
    C_k_bin = interp1(kappi,C_ki,kappa_bin,'spline');
    p_k_bin = interp1(kappi,p_ki,kappa_bin,'spline');
    bp_k_bin = interp1(kappi,bp_ki,kappa_bin,'spline');

    C_k_binl = log10(C_k_bin);
    p_k_bini = 1./p_k_bin;
    bp_k_binl = log10(bp_k_bin);

    lTp_calc3 = zeros(n_a,1);
    lTp_calc6 = zeros(n_a,1);
    lTp_calc12 = zeros(n_a,1);
    lTp_calc24 = zeros(n_a,1);
    

    for i = 1:n_a
        alpha = ray_bin(i,:);
        beta = bin_dist(i,:);
        alpha(beta==0) = [];
        vel = vel_bin(i,:);
        vel(beta==0) = [];
        beta(beta==0) = [];
        n = length(alpha);
        if (Tp_obs6(i) > 0) && (freq6(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_calc6(i) = calc_delaytime(n,alpha,beta,freq6(i),...
                aeta_binl,C_k_binl,p_k_bini,bp_k_binl,vel);
        end
        
        if (Tp_obs3(i) > 0) && (freq3(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_calc3(i) = calc_delaytime(n,alpha,beta,freq3(i),...
                aeta_binl,C_k_binl,p_k_bini,bp_k_binl,vel);
        end
        
        if (Tp_obs12(i) > 0) && (freq12(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_calc12(i) = calc_delaytime(n,alpha,beta,freq12(i),...
                aeta_binl,C_k_binl,p_k_bini,bp_k_binl,vel);
        end
        %{
        if (Tp_obs24(i) > 0) && (freq24(i) > 0) && (dist(i) >= 50) && (dist(i) <= 250)
            lTp_calc24(i) = calc_delaytime(n,alpha,beta,freq24(i),...
                aeta_binl,C_k_binl,p_k_bini,bp_k_binl,vel);
        end
        %}
    end


    Tp_error1 = sum(((Tp_obs3)-(10.^lTp_calc3)).^2);
    Tp_error2 = sum(((Tp_obs6)-(10.^lTp_calc6)).^2);
    Tp_error3 = sum(((Tp_obs12)-(10.^lTp_calc12)).^2);
    %Tp_error4 = sum(((Tp_obs24)-(10.^lTp_calc24)).^2);

    Tp_error = Tp_error2 + Tp_error1 + Tp_error3;

    
    post_prob(q) = -(Tp_error);
end

max(post_prob)
ind=find(post_prob==max(post_prob));
parm_set(ind-1,:)
parm_set(ind,:)
parm_set(ind+1,:)
parm_set(ind-9,:)
parm_set(ind+9,:)