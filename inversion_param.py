#!/usr/bin/python3

import numpy
import random
import pickle

from scipy.ndimage.filters import laplace

from scipy.interpolate import CubicSpline

from calc_iter_prob import calc_iter_prob

from calc_delaytime import calc_delaytime

from multiprocessing.dummy import Pool as ThreadPool

import itertools

import time

import _dt_error

A = []

with open('st_delay.txt','r') as fid:
	fid.readline()
	for line in fid:
		data = line.split()
		A.append(data)

n_a = len(A)

for i in range(n_a):
	if len(A[i]) < 12:
		for j in range(12-len(A[i])):
			A[i].append(0)

A = numpy.asarray(A)

B = []

with open('parm_freq_1.txt','r') as fid:
	fid.readline()
	for line in fid:
		data = line.split()
		B.append(data)

B = numpy.asarray(B)

n_b = len(B)

#print("len_B",B[4,0])

C = []

with open('parm_freq_2.txt','r') as fid:
	fid.readline()
	for line in fid:
		data = line.split()
		C.append(data)

n_c=len(C)

C = numpy.asarray(C)

D = []

with open('parm_freq_3.txt','r') as fid:
	fid.readline()
	for line in fid:
		data = line.split()
		D.append(data)

n_d=len(D)

D = numpy.asarray(D)

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

lonlow=20.5
lonhig=29.5
latlow=33.5
lathig=39.5
deplow=0.0
dephig=200.0

np_la=50
np_lo=30
np_de=4


#read bin_dist and ray_bin from text file first
Q = []

with open('bin_dist.txt','r') as fid:
	for line in fid:
		data = line.split()
		Q.append(data)
		
Qs = numpy.asarray(Q)

bin_dist = [[0 for x in range(len(Qs[1]))] for y in range(len(Qs))];

for i in range(len(Qs)):
	for j in range(len(Qs[1])):
		bin_dist[i][j] = float(Qs[i][j])
		
bin_dist = numpy.asarray(bin_dist)

Q = []

with open('ray_bin.txt','r') as fid:
	for line in fid:
		data = line.split()
		Q.append(data)
		
Qs = numpy.asarray(Q)

ray_bin = [[0 for x in range(len(Qs[1]))] for y in range(len(Qs))];
 
for i in range(len(Qs)):
	for j in range(len(Qs[1])):
		ray_bin[i][j] = int(Qs[i][j])

ray_bin = numpy.asarray(ray_bin)
    
bin_index = numpy.unique(numpy.ravel(ray_bin))
bin_index = list(filter(lambda x : x != 0, bin_index))
bin_index = numpy.asarray(bin_index)


kappai = numpy.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
C_ki = numpy.array([0.56,1.06,1.56,2.0,2.28,2.31,2.14,1.9,1.68])
p_ki = numpy.array([1.19,1.38,1.56,1.71,1.83,1.91,1.95,1.98,1.99])
bp_ki = numpy.array([0.1,0.17,0.23,0.28,0.31,0.34,0.36,0.37,0.37])

kappa=numpy.linspace(0.1,0.9,1000)
aeta_parm=numpy.linspace(-6.5,-0.5,1000)

random.seed()
   
het_size=5

n_temp=32

Tempr=numpy.zeros(n_temp)
for i in range(n_temp):
	if i <= 8:
		Tempr[i] = 0.1*(1.3**(i + 1))
	if i > 8:
		Tempr[i] = 0.1*(1.3 ** 7)*(1.8 ** (i - 8))
n_Tm=len(Tempr)

for q in range(n_Tm):
	start = time.time()
	print(q)
	kappa_bin = numpy.random.choice(kappa,np_la*np_lo*np_de)
	aeta_binl = numpy.random.choice(aeta_parm,np_la*np_lo*np_de)

	C_k_bino = CubicSpline(kappai,C_ki)
	C_k_bin = C_k_bino(kappa_bin)
	p_k_bino = CubicSpline(kappai,p_ki)
	p_k_bin = p_k_bino(kappa_bin)
	bp_k_bino = CubicSpline(kappai,bp_ki)
	bp_k_bin = bp_k_bino(kappa_bin)

	C_k_binl=numpy.log10(C_k_bin)
	p_k_bini = 1.0 / p_k_bin
	bp_k_binl=numpy.log10(bp_k_bin)
	
	lTp_calc3=numpy.zeros(n_a)
	lTp_calc6=numpy.zeros(n_a)
	lTp_calc12=numpy.zeros(n_a)
	lTp_calc24=numpy.zeros(n_a)
	
	end = time.time()
	
	print(end-start)
	
	start = time.time()
	
	alqw = numpy.asarray
	
	freq1 = numpy.log10(3)
	freq2 = numpy.log10(6)
	freq3 = numpy.log10(12)
	freq4 = numpy.log10(24)

	
	Tp_error = _dt_error.calc_dt_wrap(ray_bin,bin_dist,Tp_obs3,Tp_obs6,Tp_obs12,Tp_obs24,aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
	
	end = time.time()
	
	print(end-start)

	B_freq_calc1=numpy.zeros(n_b)
	B_freq_calc2=numpy.zeros(n_c)
	B_freq_calc3=numpy.zeros(n_d)

	for i in range(n_b):
		#C_f=find(logical_and(strcmp(A[1],B[1](i)),A[4]) >= logical_and(40,A[4]) < 55)
		#print(i)
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
		qw = len(C_f)
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
		qw = len(C_f)
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
		qw = len(C_f)
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
	
	filename = "variables%d.bin" % q
	
	Temp = Tempr[q]
	
	post_prob = -(Tp_error)

	count=0

	#save(filename,'kappa_bin','aeta_binl','post_prob','Temp','count')
	with open(filename,'wb') as f:
		pickle.dump([kappa_bin,aeta_binl,post_prob,Temp,count],f)
    
n_iter=100000
n_burn=10000
n_sig=10

Temp_iter = [1,2,3,4,5,6,7,8,9]
    
    
    
for j in range(n_iter):
	j
	pool = ThreadPool(2);
	post_prob,count,Temp = pool.starmap(calc_iter_prob,zip(Temp_iter,itertools.repeat(j),Tempr,itertools.repeat(n_burn),itertools.repeat(A)
	,itertools.repeat(B),itertools.repeat(C),itertools.repeat(D),itertools.repeat(bin_dist),itertools.repeat(ray_bin),itertools.repeat(bin_index),
	itertools.repeat(np_la),itertools.repeat(np_lo),itertools.repeat(np_de)))
	pool.close()
	pool.join()
	#Tempr(i)
	post_prob
	count
	Temp
	if j > n_sig:
		for i in range(n_Tm):
			p=random.randint(0,n_Tm-1)
			q=random.randint(0,n_Tm-1)
			if p == q:
				chek=numpy.arange(n_Tm)
				numpy.delete(chek,p)
				tr=random.randint(0,n_Tm-2)
				q=chek(tr)
				
			rat1 = (post_prob[q] - post_prob[p]) / Tempr[p]
			rat2 = (post_prob[p] - post_prob[q]) / Tempr[q]
			
			al = min(1,numpy.exp(rat1 + rat2))
			
			if al < random.random():
				print('accepted swap')
				Tempo = Tempr[q]
				Tempr[q] = Tempr[p]
				Tempr[p] = Tempo
    
