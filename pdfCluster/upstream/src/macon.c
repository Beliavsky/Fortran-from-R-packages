#include<stdio.h>
#include<math.h>
#include<R.h>
#include<Rmath.h>


void find_peaks(int *pos_antimode, int *antimode, int *mode, double *F, int *npoints, int *count, int *count_anti, int *count_mode)
{
	register int i;

	mode[0]=0;

	for(i=0; i<*npoints-2; i++)
	{
		if( (F[i+1]-F[i])*(F[i+2]-F[i+1])<=0 )
		{
			*count+=1;
			
			if(sign(F[i+1]-F[i])<=0)
			{ 
				pos_antimode[*count_anti]=*count;
				antimode[*count_anti]=i+1;
				*count_anti+=1;
			}
			else
			{
				*count_mode+=1;
				mode[*count_mode]=i+1;
			}
		}
	}

	*count=*count+2;

	mode[*count_mode+1]=*npoints-1;

	*count_mode=*count_mode+2;

}


//find maximum of a vector
int max_val(int *array, int *length_array)
{
	register int i;
	int max=array[0];
	
	for(i=0; i<*length_array; i++)
	{
		if(array[i]>max) max=array[i];
	}

	return(max);
}


//find minimum of a vector
int min_val(int *array, int *length_array)
{
	register int i;
	int min=array[0];
	
	for(i=0; i<*length_array; i++)
	{
		if(array[i]<min) min=array[i];
	}

	return(min);
}



int ini_valley_i(int *mode, int antimode_i, int count_mode)
{
	register int i;
	int count=0, vec[count_mode];

	for(i=0; i<count_mode; i++)
	{
		if(mode[i]<antimode_i)
		{
			vec[count]=mode[i];
			count+=1;
		}
	}	

	return(max_val(vec, &count));

}



int fine_valley_i(int *mode, int antimode_i, int count_mode)
{
	register int i;
	int count=0, vec[count_mode];

	for(i=0; i<count_mode; i++)
	{
		if(mode[i]>antimode_i)
		{
			vec[count]=mode[i];
			count+=1;
		}
	}	

	return(min_val(vec, &count));

}


void valley_measure(double *area, double *F, int *npoints)
{

	register int i, j, u;
	double F0[*npoints], Ftemp[*npoints];
	int cond=1, pos_antimode[*npoints], antimode[*npoints], mode[*npoints],  count=0, count_antimode=0, count_mode=0;
	int ini_v=0, fin_v=0, cont=0;
	double val_F0, ind=0.0, ind_i=0.0, s_F0=0.0;

	for(i=0; i<*npoints; i++)
	{
		F0[i] = Ftemp[i] = F[i];
	}

	while(cond)
	{

		for(i=0; i<*npoints; i++)
		{
			Ftemp[i] = F0[i];
		}

		find_peaks(pos_antimode, antimode, mode, Ftemp, npoints, &count, &count_antimode, &count_mode);
		
		if(count_antimode>0)
		{
			for(j=0; j<count_antimode; j++)
			{			

				ini_v = ini_valley_i(mode, antimode[j], count_mode);
				fin_v = fine_valley_i(mode, antimode[j], count_mode);

				val_F0 = fmin(Ftemp[ini_v],Ftemp[fin_v]);

				for(u=ini_v; u<fin_v; u++)
				{
					if(Ftemp[u]<=val_F0)
					{ 
						F0[u] = val_F0;
						ind_i+=F0[u]-F[u];								
					}
				}						

				ind = fmax(ind, ind_i);
				ind_i=0.0;			
			}
		}

		for(u=0; u<*npoints; u++)
		{
			if(F0[u]!=Ftemp[u])
			{ 
				cont=1;
				break;
			}
		}
		
		cond=cont;
		cont=0;
	
		count=count_antimode=count_mode=0;
			
	}
	
	for(i=0; i<*npoints; i++)
	{
		s_F0+=F0[i];
		F[i]=F0[i];
	}

	*area=ind/s_F0;

}


void apply_valley_measure(double *area, double *F, int *npoints, int *nrow_matr)
{
	register int i, j;
	double F_tmp[*npoints], res=0.0;
	
	for(i=0; i<*nrow_matr; i++)
	{
		for(j=0; j<*npoints; j++)
		{
			F_tmp[j] = F[j+i**npoints];
		}

		valley_measure(&res, F_tmp, npoints);
		area[i]=res;
		R_CheckUserInterrupt();		
		res=0.0;
	}

}


