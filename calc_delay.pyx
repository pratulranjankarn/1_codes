#!/usr/bin/cython

cimport numpy as np

np.import_array()

#cdef extern from "calc_delay.h":
#	double calc_delaytime(int index2,int alpha[],double beta[],int freq,double aeta_binl[],double C_k_binl[],double p_k_bini[],double bp_k_binl[]);

cdef extern from "calc_delay.h":
	double calc_delaytime(int index2,int alpha[],double beta[],int freq,double aeta_binl[],double C_k_binl[],double p_k_bini[],double bp_k_binl[]);
	double calc_delay_error(int ray_bin[15729][25],double bin_dist[15729][25],double Tp_obs3[15729],double Tp_obs6[15729],double Tp_obs12[15729],double Tp_obs24[15729],double aeta_binl[],double C_k_binl[],double p_k_bini[],double bp_k_binl[]); 

def calc_dt_wrap(np.ndarray[int, ndim=2, mode="c"] ray_bin not None,np.ndarray[double, ndim=2, mode="c"] bin_dist not None
,np.ndarray[double, ndim=1, mode="c"] Tp_obs3 not None,np.ndarray[double, ndim=1, mode="c"] Tp_obs6 not None
,np.ndarray[double, ndim=1, mode="c"] Tp_obs12 not None,np.ndarray[double, ndim=1, mode="c"] Tp_obs24 not None
,np.ndarray[double, ndim=1, mode="c"] aeta_binl not None,np.ndarray[double, ndim=1, mode="c"] C_k_binl not None
,np.ndarray[double, ndim=1, mode="c"] p_k_bini not None,np.ndarray[double, ndim=1, mode="c"] bp_k_binl not None):
	return calc_delay_error(<int(*)[25]> np.PyArray_DATA(ray_bin),<double(*)[25]> np.PyArray_DATA(bin_dist)
	,<double*> np.PyArray_DATA(Tp_obs3),<double*> np.PyArray_DATA(Tp_obs6),<double*> np.PyArray_DATA(Tp_obs12)
	,<double*> np.PyArray_DATA(Tp_obs24),<double*> np.PyArray_DATA(aeta_binl),<double*> np.PyArray_DATA(C_k_binl)
	,<double*> np.PyArray_DATA(p_k_bini),<double*> np.PyArray_DATA(bp_k_binl));
