#!/usr/bin/python3

from numpy import log10
import time

def calc_delaytime(index2,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl):
	t1 = time.time()
	if index2 == 0:
		mw=alpha[0]
		qw=beta[0]
		Mp=1 + 2*p_k_bini[mw]
		Npl=bp_k_binl[mw] + 2*p_k_bini[mw]*C_k_binl[mw] - 1.29533 + 0.18242*p_k_bini[mw] + 2*(1 - p_k_bini[mw])*aeta_binl[mw]
		lTp=Npl + Mp*log10(qw) + (2*Mp - 4)*log10(freq)
	else:
		mw=alpha[index2]
		qw=beta[index2]
		#print(mw)
		Mp=1 + 2*p_k_bini[mw]
		Npl=bp_k_binl[mw] + 2*p_k_bini[mw]*C_k_binl[mw] - 1.29533 + 0.18242*p_k_bini[mw] + 2*(1 - p_k_bini[mw])*aeta_binl[mw]
		R_dash=(calc_delaytime(index2 - 1,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl) - Npl - (2*Mp - 4)*log10(freq)) / Mp
		lTp=Npl + (2*Mp - 4)*log10(freq) + Mp*log10(10 ** R_dash + qw - beta[index2 - 1])
	
	print(t1-time.time())
	return lTp
    
