
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double calc_delaytime(int index2,int alpha[],double beta[],int freq,double aeta_binl[],double C_k_binl[],double p_k_bini[],double bp_k_binl[])
{
	int mw,sw;
	double qw, Mp, Npl, R_dash, lTp;
	if (index2 == 0)
	{
		mw=alpha[0];
		qw=beta[0];
		Mp=1 + 2*p_k_bini[mw];
		Npl=bp_k_binl[mw] + 2*p_k_bini[mw]*C_k_binl[mw] - 1.29533 + 0.18242*p_k_bini[mw] + 2*(1 - p_k_bini[mw])*aeta_binl[mw];
		lTp=Npl + Mp*log10(qw) + (2*Mp - 4)*log10(freq);
	}
	else
	{
		mw=alpha[index2];
		qw=beta[index2];
		Mp=1 + 2*p_k_bini[mw];
		Npl=bp_k_binl[mw] + 2*p_k_bini[mw]*C_k_binl[mw] - 1.29533 + 0.18242*p_k_bini[mw] + 2*(1 - p_k_bini[mw])*aeta_binl[mw];
		R_dash=(calc_delaytime(index2 - 1,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl) - Npl - (2*Mp - 4)*log10(freq)) / Mp;
		sw = index2-1;
		lTp=Npl + (2*Mp - 4)*log10(freq) + Mp*log10(pow(10,R_dash) + qw - beta[sw]);
	}
	return lTp;
}



double calc_delay_error(int ray_bin[15729][25],double bin_dist[15729][25],double Tp_obs3[15729],double Tp_obs6[15729],double Tp_obs12[15729],double Tp_obs24[15729],double aeta_binl[],double C_k_binl[],double p_k_bini[],double bp_k_binl[])
{
	int *alpha;
	double *beta;
	int i,j,k,freq;
	double Tp_error1, Tp_error2, Tp_error3, Tp_error4, Tp_error;
	double lTp_calc3[15729],lTp_calc6[15729],lTp_calc12[15729],lTp_calc24[15729];
	Tp_error1 = 0;
	Tp_error2 = 0;
	Tp_error3 = 0;
	Tp_error4 = 0;
	k=0;
	for(i=0; i<15729; i++)
	{
		for(j=0; j<25; j++)
		{
			if (bin_dist[i][j] > 0)
			{
				*alpha = ray_bin[i][j];
				*beta = bin_dist[i][j];
				alpha++;
				beta++;
				k++;
			}
		}
		freq = 3;
		if (Tp_obs3[i] > 0)
		{
			lTp_calc3[i] = calc_delaytime(k,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
		}
		freq = 6;
		if (Tp_obs6[i] > 0)
		{
			lTp_calc6[i] = calc_delaytime(k,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
		}
		freq = 12;
		if (Tp_obs12[i] > 0)
		{
			lTp_calc12[i] = calc_delaytime(k,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
		}
		freq = 24;
		if (Tp_obs24[i] > 0)
		{
			lTp_calc24[i] = calc_delaytime(k,alpha,beta,freq,aeta_binl,C_k_binl,p_k_bini,bp_k_binl);
		} 
	Tp_error1 = Tp_error1 + pow((Tp_obs3[i] - pow(10,lTp_calc3[i])),2);
	Tp_error2 = Tp_error2 + pow((Tp_obs6[i] - pow(10,lTp_calc6[i])),2);
	Tp_error3 = Tp_error3 + pow((Tp_obs12[i] - pow(10,lTp_calc12[i])),2);
	Tp_error4 = Tp_error4 + pow((Tp_obs24[i] - pow(10,lTp_calc24[i])),2);
}
Tp_error = Tp_error1 + Tp_error2 + Tp_error3 + Tp_error4;

return Tp_error;

}

