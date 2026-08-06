extern void dwt(double *Vin, int *M, int *L, double *h, double *g, 
		double *Wout, double *Vout);

extern void idwt(double *Win, double *Vin, int *M, int *L, double *h, 
		 double *g, double *Xout);

extern void modwt(double *Vin, int *N, int *j, int *L, double *ht, 
		  double *gt, double *Wout, double *Vout);

extern void imodwt(double *Win, double *Vin, int *N, int *j, int *L, 
		   double *ht, double *gt, double *Vout);
