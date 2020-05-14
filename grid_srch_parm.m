function [post_prob,count] = grid_srch_parm(kappa,aeta,A,B,C,D,bin_dist,ray_bin,bin_index,np_la,np_lo,np_de)
    %tic
    %i = 1;
    %iter = 1;
    %Temp = 1;
    %n_burn = 1;
    
    n_b = length(B{2});

    n_c = length(C{2});
     
    n_d = length(D{2});
    
    Tp_obs3 = A{7};
    Tp_obs3(isnan(Tp_obs3)) = 1;
    Tp_obs3(Tp_obs3 < 0) = 1;
    lTp_obs3 = log10(Tp_obs3);
    N_ray3 = length(Tp_obs3);

    Tp_obs6 = A{8};
    Tp_obs6(isnan(Tp_obs6)) = 1;
    Tp_obs6(Tp_obs6 < 0) = 1;
    lTp_obs6 = log10(Tp_obs6);
    N_ray6 = length(Tp_obs6);

    Tp_obs12 = A{9};
    Tp_obs12(isnan(Tp_obs12)) = 1;
    Tp_obs12(Tp_obs12 < 0) = 1;
    lTp_obs12 = log10(Tp_obs12);
    N_ray12 = length(Tp_obs12);

    Tp_obs24 = A{12};
    Tp_obs24(isnan(Tp_obs24)) = 1;
    Tp_obs24(Tp_obs24 < 0) = 1;
    lTp_obs24 = log10(Tp_obs24);
    N_ray24 = length(Tp_obs24);

    n_a = length(A{2});
    
    n_kp = length(bin_index);
    
    n_kap = 20;
    n_aet = 50;
    n_com = n_kap*n_aet;
    
    kappa = linspace(0.1,0.9,20);
    aeta = linspace(0.001,1.04,50);
    
    k =1;
    for i = 1:length(kappa)
        for j = 1: length(aeta)
            kaae(k,:) = [kappa(i),aeta(j)];
            k = k+1;
        end
    end

    
    kappai = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
    C_ki = [0.56 1.06 1.56 2.00 2.28 2.31 2.14 1.90 1.68];
    p_ki = [1.19 1.38 1.56 1.71 1.83 1.91 1.95 1.98 1.99];
    bp_ki = [0.10 0.17 0.23 0.28 0.31 0.34 0.36 0.37 0.37];
    
    

    lTp_calc3 = zeros(n_a,1);
    lTp_calc6 = zeros(n_a,1);
    lTp_calc12 = zeros(n_a,1);
    lTp_calc24 = zeros(n_a,1);


    for i = 1:n_a
      for j = 1:n_com
        alpha = ray_bin(i,:);
        alpha(alpha==0) = [];
        beta = bin_dist(i,:);
        beta(beta==0) = [];
        C_k_bin = interp1(kappai,C_ki,kaae(j,1),'spline');
        p_k_bin = interp1(kappai,p_ki,kaae(j,1),'spline');
        bp_k_bin = interp1(kappai,bp_ki,kaae(j,1),'spline');

        aeta_binl = log10(kaae(j,2));
        C_k_binl = log10(C_k_bin);
        p_k_bini = 1./p_k_bin;
        bp_k_binl = log10(bp_k_bin);
        freq = 6;
        n = length(alpha);
        if (Tp_obs6(i) > 0)
            lTp_calc6(i) = calc_delaytime(n,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
        end
        if ~isreal(lTp_calc6(i))
            %ae = aeta_binl(alpha)
            %ck = C_k_binl(alpha)
            %pk = p_k_bini(alpha)
            %bp = bp_k_binl(alpha)
            lTp_calc6(i) = 0;
        end
      end
    end

    qe1 =  find(lTp_obs6 == 0);
    lTp_calc6(qe1==0) = [];
    lTp_obs6(qe1==0) = [];
    
    qe1 =  find(lTp_obs3 == 0);
    lTp_calc3(qe1==0) = [];
    lTp_obs3(qe1==0) = [];
    
    qe1 =  find(lTp_obs12 == 0);
    lTp_calc12(qe1==0) = [];
    lTp_obs12(qe1==0) = [];
    
    qe1 =  find(lTp_obs24 == 0);
    lTp_calc24(qe1==0) = [];
    lTp_obs24(qe1==0) = [];


    Tp_error1 = sum(((10.^lTp_obs3)-(10.^lTp_calc3)).^2);
    Tp_error2 = sum(((10.^lTp_obs6)-(10.^lTp_calc6)).^2);
    Tp_error3 = sum(((10.^lTp_obs12)-(10.^lTp_calc12)).^2);
    Tp_error4 = sum(((10.^lTp_obs24)-(10.^lTp_calc24)).^2);

    Tp_error = Tp_error2;
    
    %lTp_calc6
    
    %A_f = [A{2} A{3} A{4} A{5} A{6} lTp_calc3 lTp_calc6 lTp_calc12 lTp_calc24];

    B_freq_calc1 = zeros(length(B{2}),1);
    B_freq_calc2 = zeros(length(C{2}),1);
    B_freq_calc3 = zeros(length(D{2}),1);
    
    
    for i = 1:n_b
        C_f1 = find(strcmp(A{1},B{1}(i)) & A{4} >= 40 & A{4} < 55);
        qw = length(C_f1);
        if qw > 0
            Rat3 = lTp_calc3(C_f1)-lTp_calc6(C_f1);
            Rat6 = ones(qw,1);
            Rat12 = lTp_calc12(C_f1)-lTp_calc6(C_f1);
            Rat24 = lTp_calc24(C_f1)-lTp_calc6(C_f1);
            Rat_1 = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat_1) | isnan(Rat_1) | Rat_1 == 0);
            Rat_1(mw) = [];
            freq = [ones(qw,1)*log10(3);ones(qw,1)*log10(6);ones(qw,1)*log10(12);ones(qw,1)*log10(24)];
            freq(mw) = [];
            %if ~isreal(Rat_1)
            %    save('Workspace.mat','Rat_1','-append');
            %end
            freqfit_1 = [ones(length(Rat_1),1) freq];
            parm1 = pinv(freqfit_1,0)*Rat_1;
            B_freq_calc1(i) = parm1(2);
        end
    end

    for i = 1:n_c
        C_f2 = find(strcmp(A{1},C{1}(i)) & A{4} >= 55 & A{4} < 75);
        qw = length(C_f2);
        if qw > 0
            Rat3 = lTp_calc3(C_f2)-lTp_calc6(C_f2);
            Rat6 = ones(qw,1);
            Rat12 = lTp_calc12(C_f2)-lTp_calc6(C_f2);
            Rat24 = lTp_calc24(C_f2)-lTp_calc6(C_f2);
            Rat_2 = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat_2) | isnan(Rat_2) | Rat_2 == 0);
            Rat_2(mw) = [];
            freq = [ones(qw,1)*log10(3);ones(qw,1)*log10(6);ones(qw,1)*log10(12);ones(qw,1)*log10(24)];
            freq(mw) = [];
            freqfit_2 = [ones(length(Rat_2),1) freq];
            %if ~isreal(Rat_2)
            %    save('Workspace.mat','Rat_2','-append');
            %end
            parm2 = pinv(freqfit_2,0)*Rat_2;
            B_freq_calc2(i) = parm2(2);
        end
    end
    
    for i = 1:n_d
        C_f3 = find(strcmp(A{1},D{1}(i)) & A{4} >= 75 & A{4} < 250);
        qw = length(C_f3);
        if qw > 0
            Rat3 = lTp_calc3(C_f3)-lTp_calc6(C_f3);
            Rat6 = ones(qw,1);
            Rat12 = lTp_calc12(C_f3)-lTp_calc6(C_f3);
            Rat24 = lTp_calc24(C_f3)-lTp_calc6(C_f3);
            Rat_3 = [Rat3;Rat6;Rat12;Rat24];
            mw = find(isinf(Rat_3) | isnan(Rat_3) | Rat_3 == 0);
            Rat_3(mw) = [];
            freq = [ones(qw,1)*log10(3);ones(qw,1)*log10(6);ones(qw,1)*log10(12);ones(qw,1)*log10(24)];
            freq(mw) = [];
            %if ~isreal(Rat_3)
            %    save('Workspace.mat','Rat_3','-append');
            %end
            freqfit_3 = [ones(length(Rat_3),1) freq];
            parm3 = pinv(freqfit_3,0)*Rat_3;
            B_freq_calc3(i) = parm3(2);
        end
    end

    calc1 = sum((B_freq_calc1-B{4}).^2./B{5}.^2)/n_b;
    calc2 = sum((B_freq_calc2-C{4}).^2./C{5}.^2)/n_c;
    calc3 = sum((B_freq_calc3-D{4}).^2./D{5}.^2)/n_d;

    calc_bf = w_bfreq * N_ray6 * (calc1 + calc2 + calc3);

    sq1 = 1:np_de:np_la*np_lo*np_de;
    sq2 = 2:np_de:np_la*np_lo*np_de;
    sq3 = 3:np_de:np_la*np_lo*np_de;

    Lkp1 = sum(del2(kappa_bin2(sq1)).^2);
    Lkp2 = sum(del2(kappa_bin2(sq2)).^2);
    Lkp3 = sum(del2(kappa_bin2(sq3)).^2);

    L_kp = w_kp * (Lkp1 + Lkp2 + Lkp3)/(np_la*np_lo*np_de);
    
    L_e1 = sum(del2((2./(p_k_bin(sq1)-1)).*aeta_binl(sq1)*(-0.69897)).^2);
    L_e2 = sum(del2((2./(p_k_bin(sq2)-1)).*aeta_binl(sq2)*(-0.69897)).^2);
    L_e3 = sum(del2((2./(p_k_bin(sq3)-1)).*aeta_binl(sq3)*(-0.69897)).^2);

    L_e = w_aeta * (L_e1 + L_e2 + L_e3)/(np_la*np_lo*np_de);

    L = N_ray6 * (L_e + L_kp);
    
    %calc_bf
    %Tp_error
    %L
       
    post_prob1 = -(calc_bf + L + Tp_error);
    
    paccept = ((post_prob1*Temp)-(post_prob*Tempn))/(Temp*Tempn);
    
    accep_prob = min(1,exp(paccept));
    
    pchek = rand(1);
    
    rng('shuffle');
    
    %post_prob
    
    if (accep_prob == 1 || accep_prob > rand(1))
        kappa_bin = kappa_bin2;
        aeta_bin = aeta_bin2;
        post_prob = post_prob1;
        disp('accepted');
        count = count+1
    end
    Temp = Tempn;
    save(filename,'kappa_bin','aeta_bin','post_prob','Temp','count','-append');
    
    if iter > n_burn
       if rem((iter-n_burn),10) == 0
        fid = fopen('kappa_bin.txt','a');
        fprintf(fid,'%f\t',kappa_bin);
        fclose(fid);
        fid = fopen('aeta_bin.txt','a');
        fprintf(fid,'%f\t',aeta_bin);
        fclose(fid);
       end
    end
    %toc
       %save('Workspace.mat','lTp_calc3','lTp_calc6','lTp_calc12','lTp_calc24'...
       %,'freqfit_1','freqfit_2','freqfit_3','parm1','parm2','parm3','C_f1'...
       %,'C_f2','C_f3','B_freq_calc1','B_freq_calc2','B_freq_calc3','-append');
end