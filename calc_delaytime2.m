function lTp = calc_delaytime2(beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl)
       qw = beta;
       Mp = 1+2*p_k_bini;
       Npl = bp_k_binl+2*p_k_bini*C_k_binl-1.29533+0.18242*p_k_bini+2*(1-p_k_bini)*aeta_binl;
       lTp = Npl+ Mp*log10(qw) +(2*Mp-4)*log10(freq);
end