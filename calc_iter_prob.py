#!/usr/bin/python3

import numpy
import pickle
import random

from scipy.ndimage.filters import laplace
from scipy.interpolate import CubicSpline

from calc_delaytime import calc_delaytime

    

def calc_iter_prob(i,iter,Tempn,n_burn,A,B,C,D,bin_dist,ray_bin,bin_index,np_la,np_lo,np_de):
    
	filename = "variables%d.bin" % i
    
	with open(filename,'r') as f:
		kappa_bin,aeta_binl,post_prob,Temp,count = pickle.load(f)
	
	n_b=len(B)
	n_c=len(C)
	n_d=len(D)
	n_a=len(A)
	n_kp=len(bin_index)
    
	Tp_obs3 = list(map(float,A[:,6]))
	Tp_obs3 = numpy.array(Tp_obs3)
	N_ray3=len(Tp_obs3)

	Tp_obs6 = list(map(float,A[:,7]))
	Tp_obs6 = numpy.array(Tp_obs6)
	N_ray6=len(Tp_obs6)

	Tp_obs12 = list(map(float,A[:,8]))
	Tp_obs12 = numpy.array(Tp_obs12)
	N_ray12 = len(Tp_obs12)

	Tp_obs24 = list(map(float,A[:,11]))
	Tp_obs24 = numpy.array(Tp_obs24)
	N_ray24 = len(Tp_obs24)
	
	depths = list(map(float,A[:,3]))
	depths = numpy.asarray(depths)

	w_bfreq=1
	w_kp=1
	w_aeta=1
    
	kappai = numpy.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
	C_ki = numpy.array([0.56,1.06,1.56,2.0,2.28,2.31,2.14,1.9,1.68])
	p_ki = numpy.array([1.19,1.38,1.56,1.71,1.83,1.91,1.95,1.98,1.99])
	bp_ki = numpy.array([0.1,0.17,0.23,0.28,0.31,0.34,0.36,0.37,0.37])

	random.seed()
	ind = random.randint(0,n_kp-1)
	var = bin_index[ind]
    
	kappa_bin2 = kappa_bin
	aeta_binl2 = aeta_binl
    
	random.seed()
	chus = random.randint(1,2)
	if chus == 1:
		tep=kappa_bin2[var]
		random.seed()
		kappa_bin2[var] = numpy.random.normal(kappa_bin2[var],0.05)
		while (kappa_bin2[var] > 0.9) or (kappa_bin2[var] < 0.1):
			kappa_bin2[var] = numpy.random.normal(tep,0.05)
	else:
		tep=aeta_binl2[var]
		random.seed()
		aeta_binl2[var] = numpy.random.normal(aeta_binl2[var],0.3)
		while (aeta_binl2[var] > - 0.5) or (aeta_binl2[var] < - 6.5):
			aeta_binl2[var] = numpy.random.normal(tep,0.3)


	C_k_bino = CubicSpline(kappai,C_ki)
	C_k_bin = C_k_bino(kappa_bin)
	p_k_bino = CubicSpline(kappai,p_ki)
	p_k_bin = p_k_bino(kappa_bin)
	bp_k_bino = CubicSpline(kappai,bp_ki)
	bp_k_bin = bp_k_bino(kappa_bin)

	C_k_binl = numpy.log10(C_k_bin)
	p_k_bini = 1.0 / p_k_bin
	bp_k_binl = numpy.log10(bp_k_bin)
	
	lTp_calc3=numpy.zeros(n_a)
	lTp_calc6=numpy.zeros(n_a)
	lTp_calc12=numpy.zeros(n_a)
	lTp_calc24=numpy.zeros(n_a)

	for i in range(n_a):
		alpha=ray_bin[i,:]
		alpha = list(filter(lambda x : x != 0, alpha))
		beta=bin_dist[i,:]
		beta = list(filter(lambda x : x != 0, beta))
		alpha = numpy.asarray(alpha)
		beta = numpy.asarray(beta)
		
		freq=6
		n=len(alpha)
		print(n)
		if (Tp_obs6[i] > 0):
			lTp_calc6[i]=calc_delaytime(n-1,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl)
        
		freq=3
		if (Tp_obs3[i] > 0):
			lTp_calc3[i]=calc_delaytime(n-1,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl)
            
		freq=12
		if (Tp_obs12[i] > 0):
			lTp_calc12[i]=calc_delaytime(n-1,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl)
            
		freq=24
		if (Tp_obs24[i] > 0):
			lTp_calc24[i]=calc_delaytime(n-1,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl)
            
	Tp_error1=sum(((Tp_obs3) - (10.0 ** lTp_calc3)) ** 2)
	Tp_error2=sum(((Tp_obs6) - (10.0 ** lTp_calc6)) ** 2)
	Tp_error3=sum(((Tp_obs12) - (10.0 ** lTp_calc12)) ** 2)
	Tp_error4=sum(((Tp_obs24) - (10.0 ** lTp_calc24)) ** 2)
	Tp_error=Tp_error2 + Tp_error4 + Tp_error1 + Tp_error3


	B_freq_calc1=numpy.zeros(len(B))
	B_freq_calc2=numpy.zeros(len(C))
	B_freq_calc3=numpy.zeros(len(D))

	for i in range(n_b):
		#C_f=find(logical_and(strcmp(A[1],B[1](i)),A[4]) >= logical_and(40,A[4]) < 55)
		C_f1 = (A[:,0] == B[i,0]).nonzero()
		C_f1 = numpy.asarray(C_f1)
		C_f2 = C_f1[0,:];
		C_f3 = (depths[C_f2] >= 40).nonzero()
		C_f3 = numpy.asarray(C_f3)
		C_f4 = C_f3[0,:];
		C_f5 = (depths[C_f4] < 55).nonzero()
		C_f5 = numpy.asarray(C_f5)
		C_f6 = C_f5[0,:]
		C_f = C_f2[C_f6]
		qw=len(C_f)
		if qw > 0:
			Rat3 = lTp_calc3[C_f] - lTp_calc6[C_f]
			Rat6 = numpy.zeros(qw)
			Rat12 = lTp_calc12[C_f] - lTp_calc6[C_f]
			Rat24 = lTp_calc24[C_f] - lTp_calc6[C_f]
			Rat = numpy.concatenate([Rat3,Rat6,Rat12,Rat24])
			mw = numpy.dot(1,numpy.logical_or(numpy.isinf(Rat),numpy.isnan(Rat)))
			sw = mw[mw==1]
			Rat[sw] = []
			freq = numpy.concatenate([numpy.dot(numpy.ones(qw),numpy.log10(3)),numpy.dot(numpy.ones(qw),numpy.log10(6)),
			numpy.dot(numpy.ones(qw),numpy.log10(12)),numpy.dot(numpy.ones(qw),numpy.log10(24))])
			freq[sw] = []
			freqfit = numpy.stack((numpy.ones(len(Rat)),freq))
			freqfit = numpy.transpose(freqfit)
			intrcpt,slope = numpy.linalg.lstsq(freqfit,Rat)[0]
			B_freq_calc1[i] = slope
			
	for i in range(n_c):
		#C_f=find(logical_and(strcmp(A[1],C[1](i)),A[4]) >= logical_and(55,A[4]) < 75)		
		C_f1 = (A[:,0] == C[i,0]).nonzero()
		C_f1 = numpy.asarray(C_f1)
		C_f2 = C_f1[0,:];
		C_f3 = (depths[C_f2] >= 55).nonzero()
		C_f3 = numpy.asarray(C_f3)
		C_f4 = C_f3[0,:];
		C_f5 = (depths[C_f4] < 75).nonzero()
		C_f5 = numpy.asarray(C_f5)
		C_f6 = C_f5[0,:]
		C_f = C_f2[C_f6]
		qw=len(C_f)
		if qw > 0:
			Rat3 = lTp_calc3[C_f] - lTp_calc6[C_f]
			Rat6 = numpy.zeros(qw)
			Rat12 = lTp_calc12[C_f] - lTp_calc6[C_f]
			Rat24 = lTp_calc24[C_f] - lTp_calc6[C_f]
			Rat = numpy.concatenate([Rat3,Rat6,Rat12,Rat24])
			mw = numpy.dot(1,numpy.logical_or(numpy.isinf(Rat),numpy.isnan(Rat)))
			sw = mw[mw==1]
			Rat[sw] = []
			freq = numpy.concatenate([numpy.dot(numpy.ones(qw),numpy.log10(3)),numpy.dot(numpy.ones(qw),numpy.log10(6))
			,numpy.dot(numpy.ones(qw),numpy.log10(12)),numpy.dot(numpy.ones(qw),numpy.log10(24))])
			freq[sw] = []
			freqfit = numpy.stack((numpy.ones(len(Rat)),freq))
			freqfit = numpy.transpose(freqfit)
			intrcpt,slope = numpy.linalg.lstsq(freqfit,Rat)[0]
			B_freq_calc2[i] = slope
	
	for i in range(n_d):
		#C_f=find(logical_and(strcmp(A[1],D[1](i)),A[4]) >= logical_and(75,A[4]) < 250)
		C_f1 = (A[:,0] == D[i,0]).nonzero()
		C_f1 = numpy.asarray(C_f1)
		C_f2 = C_f1[0,:];
		C_f3 = (depths[C_f2] >= 75).nonzero()
		C_f3 = numpy.asarray(C_f3)
		C_f4 = C_f3[0,:];
		C_f5 = (depths[C_f4] < 250).nonzero()
		C_f5 = numpy.asarray(C_f5)
		C_f6 = C_f5[0,:]
		C_f = C_f2[C_f6]
		qw=len(C_f)
		if qw > 0:
			Rat3 = lTp_calc3[C_f] - lTp_calc6[C_f]
			Rat6 = numpy.zeros(qw)
			Rat12 = lTp_calc12[C_f] - lTp_calc6[C_f]
			Rat24 = lTp_calc24[C_f] - lTp_calc6[C_f]
			Rat = numpy.concatenate([Rat3,Rat6,Rat12,Rat24])
			mw = numpy.dot(1,numpy.logical_or(numpy.isinf(Rat),numpy.isnan(Rat)))
			sw = mw[mw==1]
			Rat[sw] = []
			freq = numpy.concatenate([numpy.dot(numpy.ones(qw),numpy.log10(3)),numpy.dot(numpy.ones(qw),numpy.log10(6))
			,numpy.dot(numpy.ones(qw),numpy.log10(12)),numpy.dot(numpy.ones(qw),numpy.log10(24))])
			freq[sw] = []
			freqfit = numpy.stack((numpy.ones(len(Rat)),freq))
			freqfit = numpy.transpose(freqfit)
			intrcpt,slope = numpy.linalg.lstsq(freqfit,Rat)[0]
			B_freq_calc3[i] = slope
			
	B_freq_obs1 = list(map(float,B[:,3]))
	B_freq_obs1 = numpy.asarray(B_freq_obs1)
	B_freq_obs2 = list(map(float,C[:,3]))
	B_freq_obs2 = numpy.asarray(B_freq_obs2)
	B_freq_obs3 = list(map(float,D[:,3]))
	B_freq_obs3 = numpy.asarray(B_freq_obs3)
	B_freq_var1 = list(map(float,B[:,4]))
	B_freq_var1 = numpy.asarray(B_freq_var1)
	B_freq_var2 = list(map(float,C[:,4]))
	B_freq_var2 = numpy.asarray(B_freq_var2)
	B_freq_var3 = list(map(float,D[:,4]))
	B_freq_var3 = numpy.asarray(B_freq_var3)
			
	calc1 = sum(numpy.divide(numpy.power(numpy.subtract(B_freq_calc1,B_freq_obs1),2.0),B_freq_var1)) / n_b
	calc2 = sum(numpy.divide(numpy.power(numpy.subtract(B_freq_calc2,B_freq_obs2),2.0),B_freq_var2)) / n_c
	calc3 = sum(numpy.divide(numpy.power(numpy.subtract(B_freq_calc3,B_freq_obs3),2.0),B_freq_var3)) / n_d

	calc_bf = numpy.dot(numpy.dot(w_bfreq,N_ray6),(calc1 + calc2 + calc3))
    
	sq1 = numpy.arange(0,np_la*np_lo*np_de,np_de)
	sq2 = numpy.arange(1,np_la*np_lo*np_de,np_de)
	sq3 = numpy.arange(2,np_la*np_lo*np_de,np_de)
	
	Lkp1=sum(laplace(kappa_bin[sq1],mode='wrap') ** 2)
	Lkp2=sum(laplace(kappa_bin[sq2],mode='wrap') ** 2)
	Lkp3=sum(laplace(kappa_bin[sq3],mode='wrap') ** 2)
	
	L_kp = w_kp*(Lkp1 + Lkp2 + Lkp3) / (np_la*np_lo*np_de)
	
	L_e1=sum(laplace(aeta_binl[sq1],mode='wrap') ** 2)
	L_e2=sum(laplace(aeta_binl[sq2],mode='wrap') ** 2)
	L_e3=sum(laplace(aeta_binl[sq3],mode='wrap') ** 2)
	
	L_e = w_aeta*(L_e1 + L_e2 + L_e3) / (np_la*np_lo*np_de)
	
	L = N_ray6*(L_e + L_kp)
    
    
    #calc_bf
    #Tp_error
    #L
    
	post_prob1 = -(Tp_error)
	paccept=(numpy.dot(post_prob1,Temp) - numpy.dot(post_prob,Tempn)) / (Temp*Tempn)
	accep_prob=min(1,numpy.exp(paccept))
    
	random.seed()
    
    #if Tp_error < 1000:
    #    iter
    
    #post_prob
    #calc_bf
    
	if (accep_prob == 1 or accep_prob > random.random()):
		kappa_bin = kappa_bin2
		aeta_binl = aeta_binl2
		post_prob = post_prob1
		count += 1
    
	Temp = Tempn
    
	with open(filename,'wb') as f:
		pickle.dump([kappa_bin,aeta_binl,post_prob,Temp,count],f)
	
	if iter > n_burn:
		if rem((iter - n_burn),10) == 0:
			with open('kappa_binf.txt','a') as f:
				f.write('%f\\t' % (kappa_bin))
			with open('aeta_binlf.txt','a') as f:
				f.write('%f\\t' % (aeta_binl))
    
	return post_prob,count,Temp
    
