function [post_prob,count,Temp] = calc_iter_prob(inm,iter,Tempn,n_burn,...
    Tp_obs3,Tp_obs6,Tp_obs12,Bfo1,Bfv1,Bfo2,Bfv2,Bfo3,Bfv3,...
    freq3,freq6,freq12,bin_index,sq1,sq2,sq3,sq4,sq5,C_ff,C_fs,C_ft,n_l)
    %tic
    %i = 1;
    %iter = 1;
    %Temp = 1;
    %n_burn = 1;
    
    %vel = single(4.0);
    
    filename = sprintf('variables%d.mat',inm);
    
    load(filename,'kappa_bin','aeta_binl','post_prob'...
        ,'Temp','count','lTp_calc3','lTp_calc6','lTp_calc12');
    load('ray_paramfynT_15.mat','ray_bin','bin_dist','dist','vel_bin')
    
    vel_bin = 0 * single(vel_bin) + 4;
    
    n_b = single(length(Bfo1));

    n_c = single(length(Bfo2));
     
    n_d = single(length(Bfo3));
    
    N_ray6 = single(length(lTp_calc6(lTp_calc6~=0)));
    
    n_kp = single(length(bin_index));

    w_bfreq = single(3);
    w_kp = single(10);
    w_aeta = single(5);
    
    std_kp = single(0.025);
    std_ae = single(0.15);
    %std_kp2 = 0.025;
    %std_ae2 = 0.20;
    
    kappai = single([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]);
    C_ki = single([0.56 1.06 1.56 2.00 2.28 2.31 2.14 1.90 1.68 1.50]);
    p_ki = single([1.19 1.38 1.56 1.71 1.83 1.91 1.95 1.98 1.99 1.99]);
    bp_ki = single([0.10 0.17 0.23 0.28 0.31 0.34 0.36 0.37 0.37 0.38]);
    
    %rng('shuffle');
    ind = single(randi(n_kp));
    var = bin_index(ind);
    kappa_bin2 = kappa_bin;
    aeta_binl2 = aeta_binl;
    chus = randi(2);
    if chus == 1
        temp1 = kappa_bin2(var);
        temp2 = single(trandn((0.1-temp1)/std_kp,(0.9-temp1)/std_kp));
        kappa_bin2(var) = temp1 + temp2 * std_kp;
    else
        temp1 = aeta_binl2(var);
        temp2 = single(trandn((-6.5-temp1)/std_ae,(-0.5-temp1)/std_ae));
        aeta_binl2(var) = temp1 + temp2 * std_ae;
    end
    
    lTp_calcn3 = lTp_calc3;
    lTp_calcn6 = lTp_calc6;
    lTp_calcn12 = lTp_calc12;
    %lTp_calcn24 = lTp_calc24;
    
    [qerr,~] = find(ray_bin==var);
    
    qerr = single(qerr);
    
    C_k_bin = interp1(kappai,C_ki,kappa_bin2,'spline');
    p_k_bin = interp1(kappai,p_ki,kappa_bin2,'spline');
    bp_k_bin = interp1(kappai,bp_ki,kappa_bin2,'spline');
    
    C_k_binl = log10(C_k_bin);
    p_k_bini = 1./p_k_bin;
    bp_k_binl = log10(bp_k_bin);
    
    %tic
    for i = 1:length(qerr)
        qm = qerr(i);
        alpha = ray_bin(qm,:);
        beta = bin_dist(qm,:);
        vel = vel_bin(qm,:);
        vel(beta==0) = [];
        alpha(beta==0) = [];
        beta(beta==0) = [];
        
        n = single(length(alpha));
        if (dist(qm) >= 50) && (dist(qm) <= 250) && (Tp_obs6(qm) > 0)
            lTp_calcn6(qm) = calc_delaytime(n,alpha,beta,freq6,...
                aeta_binl2,C_k_binl,p_k_bini,bp_k_binl,vel);
        end
        
        if (dist(qm) >= 50) && (dist(qm) <= 250) && (Tp_obs3(qm) > 0)
            lTp_calcn3(qm) = calc_delaytime(n,alpha,beta,freq3,...
                aeta_binl2,C_k_binl,p_k_bini,bp_k_binl,vel);
        end
        
        if (dist(qm) >= 50) && (dist(qm) <= 250) && (Tp_obs12(qm) > 0) 
            lTp_calcn12(qm) = calc_delaytime(n,alpha,beta,freq12,...
                aeta_binl2,C_k_binl,p_k_bini,bp_k_binl,vel);
        end
        %{
        if (Tp_obs24(qm) > 0) && (freq24(qm) > 0) && (dist(qm) >= 50) && (dist(qm) <= 250)
            lTp_calcn24(qm) = calc_delaytime(n,alpha,beta,freq24(qm),...
                aeta_binl2,C_k_binl,p_k_bini,bp_k_binl,vel);
        end
        %}
    end
    %toc
    Tp_calc3 = 10.^lTp_calcn3;
    Tp_calc6 = 10.^lTp_calcn6;
    Tp_calc12 = 10.^lTp_calcn12;
    %Tp_calc24 = 10.^lTp_calcn24;
    
    Tp_calc3(lTp_calcn3==0) = 0;
    Tp_calc6(lTp_calcn6==0) = 0;
    Tp_calc12(lTp_calcn12==0) = 0;
    %Tp_calc24(lTp_calcn24==0) = 0;

    Tp_error1 = sum((Tp_obs3-Tp_calc3).^2)/0.04;
    Tp_error2 = sum((Tp_obs6-Tp_calc6).^2)/0.04;
    Tp_error3 = sum((Tp_obs12-Tp_calc12).^2)/0.04;
    %Tp_error4 = sum((Tp_obs24-Tp_calc24).^2)/0.04;

    Tp_errorn = Tp_error2 + Tp_error1 + Tp_error3;
    
    
    %{
    if (Tempn == 20)
        weq = ((Tp_obs24)-(10.^lTp_calc24)).^2;
        [index,werr] =  find(weq > 100);
        %index1 = find(weq==max(weq));
        %binmax = 0;
        %for qe = 1:length(index)
            %binmax = [binmax,ray_bin(qe,:)];
        %end
        %binmax = unique(binmax);
        %length(binmax)
        %length(index)
        %binmax = [];
        %binm = zeros(1,length(binmax));
        %for qr = 1:length(binmax)
        %[qerrm,qercm] = find(ray_bin==binmax(qr));
        %binm(qr) = length(qerrm);
        %end
        %hist(binm,qr)
        %kappa_bin(var)
        %((10^aeta_binl2(var))*5)^((p_k_bin(var)-1)/2)
        %length(werr)
        %Tp_obs24(index1)
    end
    %}
    
    
    %%{
    Bf_calc1 = single(zeros(n_b,1));
    Bf_calc2 = single(zeros(n_c,1));
    Bf_calc3 = single(zeros(n_d,1));
    
    
    for i = 1:n_b
        C_f1 = C_ff(i,:);
        C_f1(C_f1==0) = [];
        qw = single(length(C_f1));
        %qw
        if qw > 0
            Rat3 = lTp_calcn3(C_f1)-lTp_calcn6(C_f1);
            Rat6 = single(zeros(qw,1));
            Rat12 = lTp_calcn12(C_f1)-lTp_calcn6(C_f1);
            %Rat24 = lTp_calcn24(C_f1)-lTp_calcn6(C_f1);
            Rat_1 = [Rat3;Rat6;Rat12];
            mw = find(isinf(Rat_1) | isnan(Rat_1));
            Rat_1(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12)];
            freq(mw) = [];
            freqfit_1 = [single(ones(length(Rat_1),1)) freq];
            parm1 = pinv(freqfit_1,0)*Rat_1;
            Bf_calc1(i) = parm1(2);
        end
    end

    for i = 1:n_c
        C_f2 = C_fs(i,:);
        C_f2(C_f2==0) = [];
        qw = single(length(C_f2));
        %qw
        if qw > 0
            Rat3 = lTp_calcn3(C_f2)-lTp_calcn6(C_f2);
            Rat6 = single(zeros(qw,1));
            Rat12 = lTp_calcn12(C_f2)-lTp_calcn6(C_f2);
            %Rat24 = lTp_calcn24(C_f2)-lTp_calcn6(C_f2);
            Rat_2 = [Rat3;Rat6;Rat12];
            mw = find(isinf(Rat_2) | isnan(Rat_2));
            Rat_2(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12)];
            freq(mw) = [];
            freqfit_2 = [single(ones(length(Rat_2),1)) freq];
            parm2 = pinv(freqfit_2,0)*Rat_2;
            Bf_calc2(i) = parm2(2);
        end
    end
    
    for i = 1:n_d
        C_f3 = C_ft(i,:);
        C_f3(C_f3==0) = [];
        qw = single(length(C_f3));
        %qw
        if qw > 0
            Rat3 = lTp_calcn3(C_f3)-lTp_calcn6(C_f3);
            Rat6 = single(zeros(qw,1));
            Rat12 = lTp_calcn12(C_f3)-lTp_calcn6(C_f3);
            %Rat24 = lTp_calcn24(C_f3)-lTp_calcn6(C_f3);
            Rat_3 = [Rat3;Rat6;Rat12];
            mw = find(isinf(Rat_3) | isnan(Rat_3));
            Rat_3(mw) = [];
            freq = [single(ones(qw,1))*log10(3);single(ones(qw,1))*log10(6);single(ones(qw,1))*log10(12)];
            freq(mw) = [];
            freqfit_3 = [single(ones(length(Rat_3),1)) freq];
            parm3 = pinv(freqfit_3,0)*Rat_3;
            Bf_calc3(i) = parm3(2);
        end
    end
    
    Bf_calc1(Bfo1==0) = 0;
    Bf_calc2(Bfo2==0) = 0;
    Bf_calc3(Bfo3==0) = 0;
    Bf_calc1(Bfv1==0) = 0;
    Bf_calc2(Bfv2==0) = 0;
    Bf_calc3(Bfv3==0) = 0;
    
    qe1 = find(Bf_calc1==0);
    Bf_calc1(qe1)=[];
    Bfo1(qe1)=[];
    Bfv1(qe1)=[];
    
    n_b = length(Bfo1);
    
    qe1 = find(Bf_calc2==0);
    Bf_calc2(qe1)=[];
    Bfo2(qe1)=[];
    Bfv2(qe1)=[];
    
    n_c = length(Bfo2);
    
    qe1 = find(Bf_calc3==0);
    Bf_calc3(qe1)=[];
    Bfo3(qe1)=[];
    Bfv3(qe1)=[];
    
    n_d = length(Bfo3);

    calc1 = sum((Bf_calc1-Bfo1).^2./Bfv1)/n_b;
    calc2 = sum((Bf_calc2-Bfo2).^2./Bfv2)/n_c;
    calc3 = sum((Bf_calc3-Bfo3).^2./Bfv3)/n_d;

    calc_bfn = w_bfreq * N_ray6 * (calc1 + calc2 + calc3);
    

    Lkp1 = sum(del2(kappa_bin2(sq1)).^2);
    Lkp2 = sum(del2(kappa_bin2(sq2)).^2);
    Lkp3 = sum(del2(kappa_bin2(sq3)).^2);
    Lkp4 = sum(del2(kappa_bin2(sq4)).^2);
    Lkp5 = sum(del2(kappa_bin2(sq5)).^2);

    L_kp = w_kp * (Lkp1/n_l(1) + Lkp2/n_l(2) + Lkp3/n_l(3) + Lkp4/n_l(4) + Lkp5/n_l(5));
    
    L_e1 = sum(del2(aeta_binl2(sq1)).^2);
    L_e2 = sum(del2(aeta_binl2(sq2)).^2);
    L_e3 = sum(del2(aeta_binl2(sq3)).^2);
    L_e4 = sum(del2(aeta_binl2(sq4)).^2);
    L_e5 = sum(del2(aeta_binl2(sq5)).^2);

    L_e = w_aeta * (L_e1/n_l(1) + L_e2/n_l(2) + L_e3/n_l(3) + L_e4/n_l(4) + L_e5/n_l(5));

    Ln = N_ray6 * (L_e + L_kp);    
    %}
    
    N_ray3 = single(length(lTp_calc3(lTp_calc3~=0)));
    N_ray12 = single(length(lTp_calc12(lTp_calc12~=0)));
    %N_ray24 = single(length(lTp_calc24(lTp_calc24~=0)));
    
    %%{
    if round(Tempn)==20
        calc_bfn
        sqrt((Tp_error1+Tp_error2+Tp_error3)/(N_ray6+N_ray3+N_ray12))
        %Ln
        L_kp/w_kp
        L_e/w_aeta
        %Tp_error1
        %Tp_error2
        %Tp_error3
        %Tp_error4
        %Tp_errorn*1
        %n_b
        %n_c
        %n_d
    end
    %}
    post_prob1 = -(Ln + Tp_errorn + calc_bfn);
    
    paccept1 = ((post_prob1*Temp)-(post_prob*Tempn))/(Temp*Tempn);
    
    accep_prob1 = min(1,exp(paccept1));
    
    %rng('shuffle');
   
    
    if (accep_prob1 == 1 || accep_prob1 > rand(1))
        kappa_bin = kappa_bin2;
        aeta_binl = aeta_binl2;
        post_prob = post_prob1;
        lTp_calc3 = lTp_calcn3;
        lTp_calc6 = lTp_calcn6;
        lTp_calc12 = lTp_calcn12;
        %lTp_calc24 = lTp_calcn24;
        if iter > n_burn
            count = count+1;
        end
        %i
    %{    
    %{else
        
        kappa_bin3 = kappa_bin;
        aeta_binl3 = aeta_binl;
        if chus == 1
            rng('shuffle');
            temp2 = trandn((0.1-kappa_bin3(var))/std_kp2,(0.9-kappa_bin3(var))/std_kp2);
            kappa_bin3(var) = kappa_bin2(var) + temp2 * std_kp2;
            qmd_mdd = pdf('Normal',kappa_bin2(var),kappa_bin3(var),std_kp);
            qmd_m = pdf('Normal',kappa_bin2(var),kappa_bin(var),std_kp);
        %disp(var);
        else
            rng('shuffle');
            temp2 = trandn((-6.5-aeta_binl3(var))/std_ae2,(-0.5-aeta_binl3(var))/std_ae2);
            aeta_binl3(var) = aeta_binl3(var) + temp2 * std_ae2;
            qmd_mdd = pdf('Normal',aeta_binl2(var),aeta_binl3(var),std_ae);
            qmd_m = pdf('Normal',aeta_binl2(var),aeta_binl(var),std_ae);
        %disp(aeta_binl(var));
        end
        C_k_bin = interp1(kappai,C_ki,kappa_bin3,'spline');
        p_k_bin = interp1(kappai,p_ki,kappa_bin3,'spline');
        bp_k_bin = interp1(kappai,bp_ki,kappa_bin3,'spline');
    
        C_k_binl = log10(C_k_bin);
        p_k_bini = 1./p_k_bin;
        bp_k_binl = log10(bp_k_bin);

        lTp_calc3 = zeros(n_a,1);
        lTp_calc6 = zeros(n_a,1);
        lTp_calc12 = zeros(n_a,1);
        lTp_calc24 = zeros(n_a,1);
    

    %tic
        for i = 1:n_a
            alpha = ray_bin(i,:);
            %alpha(alpha==0) = [];
            beta = bin_dist(i,:);
            alpha(beta==0) = [];
            beta(beta==0) = [];
        
            n = length(alpha);
            if (Tp_obs6(i) > 0) && (freq6(i) > 0)
                lTp_calc6(i) = calc_delaytime(n,alpha,beta,freq6(i),aeta_binl3,C_k_binl,p_k_bini,bp_k_binl);
            end
        
            if (Tp_obs3(i) > 0) && (freq3(i) > 0)
                lTp_calc3(i) = calc_delaytime(n,alpha,beta,freq3(i),aeta_binl3,C_k_binl,p_k_bini,bp_k_binl);
            end
        
            if (Tp_obs12(i) > 0) && (freq12(i) > 0)
                lTp_calc12(i) = calc_delaytime(n,alpha,beta,freq12(i),aeta_binl3,C_k_binl,p_k_bini,bp_k_binl);
            end
        
            if (Tp_obs24(i) > 0) && (freq24(i) > 0)
                lTp_calc24(i) = calc_delaytime(n,alpha,beta,freq24(i),aeta_binl3,C_k_binl,p_k_bini,bp_k_binl);
            end
        
        end
    %toc

        Tp_error1 = sum(((Tp_obs3)-(10.^lTp_calc3)).^2);
        Tp_error2 = sum(((Tp_obs6)-(10.^lTp_calc6)).^2);
        Tp_error3 = sum(((Tp_obs12)-(10.^lTp_calc12)).^2);
        Tp_error4 = sum(((Tp_obs24)-(10.^lTp_calc24)).^2);

        Tp_error = Tp_error2 + Tp_error1 + Tp_error3 + Tp_error4;
        
        Bf_calc1 = zeros(length(B{2}),1);
        Bf_calc2 = zeros(length(C{2}),1);
        Bf_calc3 = zeros(length(D{2}),1);
    
    
        for i = 1:n_b
            C_f1 = find(strcmp(A{1},B{1}(i)) & A{4} >= 40 & A{4} < 55);
            qw = length(C_f1);
            if qw > 0
                Rat3 = lTp_calc3(C_f1)-lTp_calc6(C_f1);
                Rat6 = zeros(qw,1);
                Rat12 = lTp_calc12(C_f1)-lTp_calc6(C_f1);
                Rat24 = lTp_calc24(C_f1)-lTp_calc6(C_f1);
                Rat_1 = [Rat3;Rat6;Rat12;Rat24];
                mw = find(isinf(Rat_1) | isnan(Rat_1));
                Rat_1(mw) = [];
                freq = [ones(qw,1)*log10(3);ones(qw,1)*log10(6);ones(qw,1)*log10(12);ones(qw,1)*log10(24)];
                freq(mw) = [];
                %if ~isreal(Rat_1)
                %    save('Workspace.mat','Rat_1','-append');
                %end
                freqfit_1 = [ones(length(Rat_1),1) freq];
                parm1 = pinv(freqfit_1,0)*Rat_1;
                Bf_calc1(i) = parm1(2);
            end
        end

        for i = 1:n_c
            C_f2 = find(strcmp(A{1},C{1}(i)) & A{4} >= 55 & A{4} < 75);
            qw = length(C_f2);
            if qw > 0
                Rat3 = lTp_calc3(C_f2)-lTp_calc6(C_f2);
                Rat6 = zeros(qw,1);
                Rat12 = lTp_calc12(C_f2)-lTp_calc6(C_f2);
                Rat24 = lTp_calc24(C_f2)-lTp_calc6(C_f2);
                Rat_2 = [Rat3;Rat6;Rat12;Rat24];
                mw = find(isinf(Rat_2) | isnan(Rat_2));
                Rat_2(mw) = [];
                freq = [ones(qw,1)*log10(3);ones(qw,1)*log10(6);ones(qw,1)*log10(12);ones(qw,1)*log10(24)];
                freq(mw) = [];
                freqfit_2 = [ones(length(Rat_2),1) freq];
                %if ~isreal(Rat_2)
                %    save('Workspace.mat','Rat_2','-append');
                %end
                parm2 = pinv(freqfit_2,0)*Rat_2;
                Bf_calc2(i) = parm2(2);
            end
        end
    
        for i = 1:n_d
            C_f3 = find(strcmp(A{1},D{1}(i)) & A{4} >= 75 & A{4} < 250);
            qw = length(C_f3);
            if qw > 0
                Rat3 = lTp_calc3(C_f3)-lTp_calc6(C_f3);
                Rat6 = zeros(qw,1);
                Rat12 = lTp_calc12(C_f3)-lTp_calc6(C_f3);
                Rat24 = lTp_calc24(C_f3)-lTp_calc6(C_f3);
                Rat_3 = [Rat3;Rat6;Rat12;Rat24];
                mw = find(isinf(Rat_3) | isnan(Rat_3));
                Rat_3(mw) = [];
                freq = [ones(qw,1)*log10(3);ones(qw,1)*log10(6);ones(qw,1)*log10(12);ones(qw,1)*log10(24)];
                freq(mw) = [];
                %if ~isreal(Rat_3)
                %    save('Workspace.mat','Rat_3','-append');
                %end
                freqfit_3 = [ones(length(Rat_3),1) freq];
                parm3 = pinv(freqfit_3,0)*Rat_3;
                Bf_calc3(i) = parm3(2);
            end
        end
    
    
    
        Bfo1 = B{4};
        Bfv1 = B{5};
        Bfo2 = C{4};
        Bfv2 = C{5};
        Bfo3 = D{4};
        Bfv3 = D{5};
    
    
        qe1 = find(Bf_calc1==0);
        Bf_calc1(qe1)=[];
        Bfo1(qe1)=[];
        Bfv1(qe1)=[];
    
        n_b = length(Bfo1);
    
        qe1 = find(Bf_calc2==0);
        Bf_calc2(qe1)=[];
        Bfo2(qe1)=[];
        Bfv2(qe1)=[];
    
        n_c = length(Bfo2);
    
        qe1 = find(Bf_calc3==0);
        Bf_calc3(qe1)=[];
        Bfo3(qe1)=[];
        Bfv3(qe1)=[];
    
        n_d = length(Bfo3);
    

        calc1 = sum((Bf_calc1-Bfo1).^2./Bfv1)/n_b;
        calc2 = sum((Bf_calc2-Bfo2).^2./Bfv2)/n_c;
        calc3 = sum((Bf_calc3-Bfo3).^2./Bfv3)/n_d;

        calc_bf = w_bfreq * N_ray6 * (calc1 + calc2 + calc3);

        Lkp1 = sum(del2(kappa_bin3(sq1)).^2);
        Lkp2 = sum(del2(kappa_bin3(sq2)).^2);
        Lkp3 = sum(del2(kappa_bin3(sq3)).^2);
        Lkp4 = sum(del2(kappa_bin3(sq4)).^2);

        L_kp = w_kp * (Lkp1 + Lkp2 + Lkp3 + Lkp4)/(length(bin_index));
    
        L_e1 = sum(del2(aeta_binl3(sq1)).^2);
        L_e2 = sum(del2(aeta_binl3(sq2)).^2);
        L_e3 = sum(del2(aeta_binl3(sq3)).^2);
        L_e4 = sum(del2(aeta_binl3(sq4)).^2);

        L_e = w_aeta * (L_e1 + L_e2 + L_e3 + L_e4)/(length(bin_index));

        L = N_ray6 * (L_e + L_kp);
        
        post_prob2 = -(Tp_error + L + calc_bf);
    
        paccept2 = ((post_prob2*Temp)-(post_prob*Tempn))/(Temp*Tempn);
    
        accep_prob2 = min(1,exp(paccept2));
        
        paccept3 = (post_prob1-post_prob2);
        
        accep_prob3 = min(1,exp(paccept3));
        
        accep_prob4 = min(1,(exp(paccept2)*(qmd_mdd/qmd_m)...
            *(1-accep_prob3)/(1-accep_prob1)));
        
        if (accep_prob4 == 1 || accep_prob4 > rand(1))
            kappa_bin = kappa_bin3;
            aeta_binl = aeta_binl3;
            post_prob = post_prob2;
            %disp('accepted');
            count = count+1;
            if (Tempn == 1024)
                disp('accepted');
            end
        end
        %}
    end
    Temp = Tempn;
    
    
    save(filename,'kappa_bin','aeta_binl','post_prob'...
        ,'Temp','count','lTp_calc3','lTp_calc6','lTp_calc12','-append');
    
    clearvars -except post_prob count Temp ;
    
 
    %toc
end
    