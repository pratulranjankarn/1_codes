function lTp = calc_delaytime(index2,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl,vel) %%codegen
    %tic
    if index2 == 1
       mw = alpha(1);
       qw = beta(1);
       Mp = 1+2*p_k_bini(mw);
       Npl = bp_k_binl(mw)+2*p_k_bini(mw)*C_k_binl(mw)-1.29533+0.78448*p_k_bini(mw)...
           +2*(1-p_k_bini(mw))*aeta_binl(mw)-4*(1-4*p_k_bini(mw))+vel(index2)*(1-4*p_k_bini(mw));
       lTp = Npl+ Mp*log10(qw) +(2*Mp-4)*log10(freq);
    else
       mw = alpha(index2);
       qw = beta(index2);
       Mp = 1+2*p_k_bini(mw);
       Npl = bp_k_binl(mw)+2*p_k_bini(mw)*C_k_binl(mw)-1.29533+0.78448*p_k_bini(mw)...
           +2*(1-p_k_bini(mw))*aeta_binl(mw)-4*(1-4*p_k_bini(mw))+vel(index2)*(1-4*p_k_bini(mw));
       R_dash = (calc_delaytime(index2-1,alpha,beta,freq,aeta_binl,C_k_binl...
           ,p_k_bini,bp_k_binl,vel)-Npl-(2*Mp-4)*log10(freq))/Mp;
       lTp = Npl + (2*Mp-4)*log10(freq)+Mp*log10(10^R_dash+qw);
    end
    %toc
end
