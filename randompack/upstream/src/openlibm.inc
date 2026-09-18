// -*- C -*-
// The functions in this file are derived from openlibm
// (https://github.com/JuliaMath/openlibm),
// which in turn is based on the fdlibm / FreeBSD msun math library originally developed
// by Sun Microsystems.
// openlibm is licensed under the BSD 2-Clause license.

#include <float.h>
#include <stdint.h>
#include <string.h>

static inline uint64_t rp_u64_from_f64(double x) {
  uint64_t u;
  memcpy(&u, &x, sizeof(u));
  return u;
}

static inline double rp_f64_from_u64(uint64_t u) {
  double x;
  memcpy(&x, &u, sizeof(x));
  return x;
}

static inline uint32_t rp_u32_from_f32(float x) {
  uint32_t u;
  memcpy(&u, &x, sizeof(u));
  return u;
}

static inline float rp_f32_from_u32(uint32_t u) {
  float x;
  memcpy(&x, &u, sizeof(x));
  return x;
}

#define EXTRACT_WORDS(ix0, ix1, d) do {         \
   uint64_t u__ = rp_u64_from_f64((d));         \
   (ix0) = (uint32_t)(u__ >> 32);               \
   (ix1) = (uint32_t)(u__ & 0xffffffffu);       \
 } while (0)

#define SET_HIGH_WORD(d, v) do {                \
   uint64_t u__ = rp_u64_from_f64((d));         \
   u__ = ((uint64_t)(uint32_t)(v) << 32) |      \
     (u__ & 0xffffffffu);                       \
   (d) = rp_f64_from_u64(u__);                  \
 } while (0)

#define GET_LOW_WORD(i, d) do {                    \
   uint64_t u__ = rp_u64_from_f64((d));            \
   (i) = (int32_t)(uint32_t)(u__ & 0xffffffffu);   \
 } while (0)

#define GET_HIGH_WORD(i, d) do {                \
   uint64_t u__ = rp_u64_from_f64((d));         \
   (i) = (int32_t)(uint32_t)(u__ >> 32);        \
 } while (0)

#define INSERT_WORDS(d, ix0, ix1) do {                \
   uint64_t u__ = ((uint64_t)(uint32_t)(ix0) << 32) | \
     (uint64_t)(uint32_t)(ix1);                       \
   (d) = rp_f64_from_u64(u__);                        \
 } while (0)

#define GET_FLOAT_WORD(i, d) do {               \
   uint32_t u__ = rp_u32_from_f32((d));         \
   (i) = (int32_t)u__;                          \
 } while (0)

#define GET_FLOAT_WORD_U(i, d) do {             \
   (i) = rp_u32_from_f32((d));                  \
 } while (0)

#define SET_FLOAT_WORD(d, i) do {               \
   (d) = rp_f32_from_u32((uint32_t)(i));        \
 } while (0)

static const double
ln2_hi  =  6.93147180369123816490e-01,	/* 3fe62e42 fee00000 */
  ln2_lo  =  1.90821492927058770002e-10,	/* 3dea39ef 35793c76 */
  two54   =  1.80143985094819840000e+16,  /* 43500000 00000000 */
  Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
  Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
  Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
  Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
  Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
  Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
  Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */
static const double
Lp1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
  Lp2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
  Lp3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
  Lp4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
  Lp5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
  Lp6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
  Lp7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */

static const double zero   =  0.0;

ALWAYS_INLINE double openlibm_log(double x)
{
  double hfsq,f,s,z,R,w,t1,t2,dk;
  int32_t k,hx,i,j;
  uint32_t lx;
  EXTRACT_WORDS(hx,lx,x);
  k=0;
  if (hx < 0x00100000) {			/* x < 2**-1022  */
    if (((hx&0x7fffffff)|lx)==0) 
		return -two54/zero;		/* log(+-0)=-inf */
    if (hx<0) return NAN;	/* log(-#) = NaN */
    k -= 54; x *= two54; /* subnormal number, scale up x */
    GET_HIGH_WORD(hx,x);
  } 
  if (hx >= 0x7ff00000) return x+x;
  k += (hx>>20)-1023;
  hx &= 0x000fffff;
  i = (hx+0x95f64)&0x100000;
  SET_HIGH_WORD(x,hx|(i^0x3ff00000));	/* normalize x or x/2 */
  k += (i>>20);
  f = x-1.0;
  if((0x000fffff&(2+hx))<3) {	/* -2**-20 <= f < 2**-20 */
    if(f==zero) {
		if(k==0) {
        return zero;
		} else {
        dk=(double)k;
        return dk*ln2_hi+dk*ln2_lo;
		}
    }
    R = f*f*(0.5-0.33333333333333333*f);
    if(k==0) return f-R; else {dk=(double)k;
      return dk*ln2_hi-((R-dk*ln2_lo)-f);}
  }
  s = f/(2.0+f); 
  dk = (double)k;
  z = s*s;
  i = hx-0x6147a;
  w = z*z;
  j = 0x6b851-hx;
  t1= w*(Lg2+w*(Lg4+w*Lg6)); 
  t2= z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7))); 
  i |= j;
  R = t2+t1;
  if(i>0) {
    hfsq=0.5*f*f;
    if(k==0) return f-(hfsq-s*(hfsq+R)); else
      return dk*ln2_hi-((hfsq-(s*(hfsq+R)+dk*ln2_lo))-f);
  } else {
    if(k==0) return f-s*(f-R); else
      return dk*ln2_hi-((s*(f-R)-dk*ln2_lo)-f);
  }
}

ALWAYS_INLINE double openlibm_log1p(double x)
{
  double hfsq,f,c,s,z,R,u;
  int32_t k,hx,hu,ax;
  GET_HIGH_WORD(hx,x);
  ax = hx&0x7fffffff;
  k = 1;
  if (hx < 0x3FDA827A) {			/* 1+x < sqrt(2)+ */
    if(ax>=0x3ff00000) {		/* x <= -1.0 */
		if(x==-1.0) return -two54/zero; /* log1p(-1)=+inf */
		else return (x-x)/(x-x);	/* log1p(x<-1)=NaN */
    }
    if(ax<0x3e200000) {			/* |x| < 2**-29 */
		if(two54+x>zero			/* raise inexact */
         &&ax<0x3c900000) 		/* |x| < 2**-54 */
        return x;
		else
        return x - x*x*0.5;
    }
    if(hx>0||hx<=((int32_t)0xbfd2bec4)) {
		k=0;f=x;hu=1;}		/* sqrt(2)/2- <= 1+x < sqrt(2)+ */
  }
  if (hx >= 0x7ff00000) return x+x;
  if(k!=0) {
    if(hx<0x43400000) {
      u = 1.0 + x; //STRICT_ASSIGN(double,u,1.0+x);
      GET_HIGH_WORD(hu,u);
      k  = (hu>>20)-1023;
      c  = (k>0)? 1.0-(u-x):x-(u-1.0);/* correction term */
		c /= u;
    } else {
		u  = x;
		GET_HIGH_WORD(hu,u);
      k  = (hu>>20)-1023;
		c  = 0;
    }
    hu &= 0x000fffff;
    if(hu<0x6a09e) {			/* u ~< sqrt(2) */
      SET_HIGH_WORD(u,hu|0x3ff00000);	/* normalize u */
    } else {
      k += 1;
		SET_HIGH_WORD(u,hu|0x3fe00000);	/* normalize u/2 */
      hu = (0x00100000-hu)>>2;
    }
    f = u-1.0;
  }
  hfsq=0.5*f*f;
  if(hu==0) {	/* |f| < 2**-20 */
    if(f==zero) {
		if(k==0) {
        return zero;
		} else {
        c += k*ln2_lo;
        return k*ln2_hi+c;
		}
    }
    R = hfsq*(1.0-0.66666666666666666*f);
    if(k==0) return f-R; else
      return k*ln2_hi-((R-(k*ln2_lo+c))-f);
  }
  s = f/(2.0+f);
  z = s*s;
  R = z*(Lp1+z*(Lp2+z*(Lp3+z*(Lp4+z*(Lp5+z*(Lp6+z*Lp7))))));
  if(k==0) return f-(hfsq-s*(hfsq+R)); else
    return k*ln2_hi-((hfsq-(s*(hfsq+R)+(k*ln2_lo+c)))-f);
}

static const double
one	= 1.0,
  halF[2]	= {0.5,-0.5,},
  huge	= 1.0e+300,
  o_threshold=  7.09782712893383973096e+02,  /* 0x40862E42, 0xFEFA39EF */
  u_threshold= -7.45133219101941108420e+02,  /* 0xc0874910, 0xD52D3051 */
  ln2HI[2]   ={ 6.93147180369123816490e-01,  /* 0x3fe62e42, 0xfee00000 */
  -6.93147180369123816490e-01,},/* 0xbfe62e42, 0xfee00000 */
  ln2LO[2]   ={ 1.90821492927058770002e-10,  /* 0x3dea39ef, 0x35793c76 */
  -1.90821492927058770002e-10,},/* 0xbdea39ef, 0x35793c76 */
  invln2 =  1.44269504088896338700e+00, /* 0x3ff71547, 0x652b82fe */
  P1   =  1.66666666666666019037e-01, /* 0x3FC55555, 0x5555553E */
  P2   = -2.77777777770155933842e-03, /* 0xBF66C16C, 0x16BEBD93 */
  P3   =  6.61375632143793436117e-05, /* 0x3F11566A, 0xAF25DE2C */
  P4   = -1.65339022054652515390e-06, /* 0xBEBBBD41, 0xC5D26BF1 */
  P5   =  4.13813679705723846039e-08; /* 0x3E663769, 0x72BEA4D0 */

static volatile double
twom1000= 9.33263618503218878990e-302;     /* 2**-1000=0x01700000,0*/

ALWAYS_INLINE double openlibm_exp(double x)
{
  double y,hi=0.0,lo=0.0,c,t,twopk;
  int32_t k=0,xsb;
  uint32_t hx;
  GET_HIGH_WORD(hx,x);
  xsb = (hx>>31)&1;		/* sign bit of x */
  hx &= 0x7fffffff;		/* high word of |x| */
  /* filter out non-finite argument */
  if(hx >= 0x40862E42) {			/* if |x|>=709.78... */
    if(hx>=0x7ff00000) {
      uint32_t lx;
		GET_LOW_WORD(lx,x);
		if(((hx&0xfffff)|lx)!=0)
        return x+x; 		/* NaN */
		else return (xsb==0)? x:0.0;	/* exp(+-inf)={inf,0} */
    }
    if(x > o_threshold) return huge*huge; /* overflow */
    if(x < u_threshold) return twom1000*twom1000; /* underflow */
  }
  if (x == 1.0)
    return 2.718281828459045235360;
  /* argument reduction */
  if(hx > 0x3fd62e42) {		/* if  |x| > 0.5 ln2 */ 
    if(hx < 0x3FF0A2B2) {	/* and |x| < 1.5 ln2 */
		hi = x-ln2HI[xsb]; lo=ln2LO[xsb]; k = 1-xsb-xsb;
    } else {
		k  = (int)(invln2*x+halF[xsb]);
		t  = k;
		hi = x - t*ln2HI[0];	/* t*ln2HI is exact here */
		lo = t*ln2LO[0];
    }
    x = hi - lo;
  } 
  else if(hx < 0x3e300000)  {	/* when |x|<2**-28 */
    if(huge+x>one) return one+x;/* trigger inexact */
  }
  else k = 0;
  /* x is now in primary range */
  t  = x*x;
  if(k >= -1021)
    INSERT_WORDS(twopk,0x3ff00000+(k<<20), 0);
  else
    INSERT_WORDS(twopk,0x3ff00000+((k+1000)<<20), 0);
  c  = x - t*(P1+t*(P2+t*(P3+t*(P4+t*P5))));
  if(k==0) 	return one-((x*c)/(c-2.0)-x); 
  else 		y = one-((lo-(x*c)/(2.0-c))-hi);
  if(k >= -1021) {
    if (k==1024) return y*2.0*0x1p1023;
    return y*twopk;
  } else {
    return y*twopk*twom1000;
  }
}

static const float
ln2_hi_f =   6.9313812256e-01f,	/* 0x3f317180 */
  ln2_lo_f =   9.0580006145e-06f,	/* 0x3717f7d1 */
  two25_f =    3.355443200e+07f,	/* 0x4c000000 */
  Lg1_f =      0xaaaaaa.0p-24f,	/* 0.66666662693 */
  Lg2_f =      0xccce13.0p-25f,	/* 0.40000972152 */
  Lg3_f =      0x91e9ee.0p-25f,	/* 0.28498786688 */
  Lg4_f =      0xf89e26.0p-26f;	/* 0.24279078841 */

static const float zero_f = 0.0f;

ALWAYS_INLINE float openlibm_logf(float x)
{
  float hfsq,f,s,z,R,w,t1,t2,dk;
  int32_t k,ix,i,j;
  GET_FLOAT_WORD(ix,x);
  k=0;
  if (ix < 0x00800000) {			/* x < 2**-126  */
    if ((ix&0x7fffffff)==0)
		return -two25_f/zero_f;		/* log(+-0)=-inf */
    if (ix<0) return NAN;	/* log(-#) = NaN */
    k -= 25; x *= two25_f; /* subnormal number, scale up x */
    GET_FLOAT_WORD(ix,x);
  }
  if (ix >= 0x7f800000) return x+x;
  k += (ix>>23)-127;
  ix &= 0x007fffff;
  i = (ix+(0x95f64<<3))&0x800000;
  SET_FLOAT_WORD(x,ix|(i^0x3f800000));	/* normalize x or x/2 */
  k += (i>>23);
  f = x-(float)1.0;
  if((0x007fffff&(0x8000+ix))<0xc000) {	/* -2**-9 <= f < 2**-9 */
    if(f==zero_f) {
		if(k==0) {
        return zero_f;
		} else {
        dk=(float)k;
        return dk*ln2_hi_f+dk*ln2_lo_f;
		}
    }
    R = f*f*((float)0.5-(float)0.33333333333333333*f);
    if(k==0) return f-R; else {dk=(float)k;
      return dk*ln2_hi_f-((R-dk*ln2_lo_f)-f);}
  }
  s = f/((float)2.0+f);
  dk = (float)k;
  z = s*s;
  i = ix-(0x6147a<<3);
  w = z*z;
  j = (0x6b851<<3)-ix;
  t1= w*(Lg2_f+w*Lg4_f);
  t2= z*(Lg1_f+w*Lg3_f);
  i |= j;
  R = t2+t1;
  if(i>0) {
    hfsq=(float)0.5*f*f;
    if(k==0) return f-(hfsq-s*(hfsq+R)); else
      return dk*ln2_hi_f-((hfsq-(s*(hfsq+R)+dk*ln2_lo_f))-f);
  } else {
    if(k==0) return f-s*(f-R); else
      return dk*ln2_hi_f-((s*(f-R)-dk*ln2_lo_f)-f);
  }
}

static const float
ln2_hi1_f =   6.9313812256e-01f,	/* 0x3f317180 */
  ln2_lo1_f =   9.0580006145e-06f,	/* 0x3717f7d1 */
  two25_1f =    3.355443200e+07f,	/* 0x4c000000 */
  Lp1_f = 6.6666668653e-01f,	/* 3F2AAAAB */
  Lp2_f = 4.0000000596e-01f,	/* 3ECCCCCD */
  Lp3_f = 2.8571429849e-01f, /* 3E924925 */
  Lp4_f = 2.2222198546e-01f, /* 3E638E29 */
  Lp5_f = 1.8183572590e-01f, /* 3E3A3325 */
  Lp6_f = 1.5313838422e-01f, /* 3E1CD04F */
  Lp7_f = 1.4798198640e-01f; /* 3E178897 */

ALWAYS_INLINE float openlibm_log1pf(float x)
{
  float hfsq,f,c,s,z,R,u;
  int32_t k,hx,hu,ax;
  GET_FLOAT_WORD(hx,x);
  ax = hx&0x7fffffff;
  k = 1;
  if (hx < 0x3ed413d0) {			/* 1+x < sqrt(2)+  */
    if(ax>=0x3f800000) {		/* x <= -1.0 */
		if(x==(float)-1.0) return -two25_1f/zero_f; /* log1p(-1)=+inf */
		else return (x-x)/(x-x);	/* log1p(x<-1)=NaN */
    }
    if(ax<0x38000000) {			/* |x| < 2**-15 */
		if(two25_1f+x>zero_f			/* raise inexact */
         &&ax<0x33800000) 		/* |x| < 2**-24 */
        return x;
		else
        return x - x*x*(float)0.5;
    }
    if(hx>0||hx<=((int32_t)0xbe95f619)) {
		k=0;f=x;hu=1;}		/* sqrt(2)/2- <= 1+x < sqrt(2)+ */
  }
  if (hx >= 0x7f800000) return x+x;
  if(k!=0) {
    if(hx<0x5a000000) {
		u = (float)1.0+x;
		GET_FLOAT_WORD(hu,u);
      k  = (hu>>23)-127;
		/* correction term */
      c  = (k>0)? (float)1.0-(u-x):x-(u-(float)1.0);
		c /= u;
    } else {
		u  = x;
		GET_FLOAT_WORD(hu,u);
      k  = (hu>>23)-127;
		c  = 0;
    }
    hu &= 0x007fffff;
    if(hu<0x3504f4) {			/* u < sqrt(2) */
      SET_FLOAT_WORD(u,hu|0x3f800000);/* normalize u */
    } else {
      k += 1;
		SET_FLOAT_WORD(u,hu|0x3f000000);	/* normalize u/2 */
      hu = (0x00800000-hu)>>2;
    }
    f = u-(float)1.0;
  }
  hfsq=(float)0.5*f*f;
  if(hu==0) {	/* |f| < 2**-20 */
    if(f==zero_f) {
		if(k==0) {
        return zero_f;
		} else {
        c += k*ln2_lo1_f;
        return k*ln2_hi1_f+c;
		}
    }
    R = hfsq*((float)1.0-(float)0.66666666666666666*f);
    if(k==0) return f-R; else
      return k*ln2_hi1_f-((R-(k*ln2_lo1_f+c))-f);
  }
  s = f/((float)2.0+f);
  z = s*s;
  R = z*(Lp1_f+z*(Lp2_f+z*(Lp3_f+z*(Lp4_f+z*(Lp5_f+z*(Lp6_f+z*Lp7_f))))));
  if(k==0) return f-(hfsq-s*(hfsq+R)); else
    return k*ln2_hi1_f-((hfsq-(s*(hfsq+R)+(k*ln2_lo1_f+c)))-f);
}

static const float
one_f	= 1.0f,
  halF_f[2]	= {0.5f,-0.5f,},
  huge_f	= 1.0e+30f,
  o_threshold_f=  8.8721679688e+01f,  /* 0x42b17180 */
  u_threshold_f= -1.0397208405e+02f,  /* 0xc2cff1b5 */
  ln2HI_f[2]   ={ 6.9314575195e-01f,		/* 0x3f317200 */
  -6.9314575195e-01f,},	/* 0xbf317200 */
  ln2LO_f[2]   ={ 1.4286067653e-06f,  	/* 0x35bfbe8e */
  -1.4286067653e-06f,},	/* 0xb5bfbe8e */
  invln2_f =  1.4426950216e+00f, 		/* 0x3fb8aa3b */
  P1_f =  1.6666625440e-1f,		/*  0xaaaa8f.0p-26 */
  P2_f = -2.7667332906e-3f;		/* -0xb55215.0p-32 */

static volatile float twom100_f = 7.8886090522e-31f;      /* 2**-100=0x0d800000 */

ALWAYS_INLINE float openlibm_expf(float x)
{
  float y,hi=0.0f,lo=0.0f,c,t,twopk;
  int32_t k=0,xsb;
  uint32_t hx;
  GET_FLOAT_WORD_U(hx,x);
  xsb = (hx>>31)&1;		/* sign bit of x */
  hx &= 0x7fffffff;		/* high word of |x| */
  /* filter out non-finite argument */
  if(hx >= 0x42b17218) {			/* if |x|>=88.721... */
    if(hx>0x7f800000)
      return x+x;	 		/* NaN */
    if(hx==0x7f800000)
		return (xsb==0)? x:0.0f;		/* exp(+-inf)={inf,0} */
    if(x > o_threshold_f) return huge_f*huge_f; /* overflow */
    if(x < u_threshold_f) return twom100_f*twom100_f; /* underflow */
  }
  /* argument reduction */
	if(hx > 0x3eb17218) {		/* if  |x| > 0.5 ln2 */
	  if(hx < 0x3F851592) {	/* and |x| < 1.5 ln2 */
			hi = x-ln2HI_f[xsb]; lo=ln2LO_f[xsb]; k = 1-xsb-xsb;
	  } else {
			k  = (int32_t)(invln2_f*x+halF_f[xsb]);
			t  = (float)k;
			hi = x - t*ln2HI_f[0];	/* t*ln2HI is exact here */
			lo = t*ln2LO_f[0];
	  }
    x = hi - lo;
  }
  else if(hx < 0x39000000)  {	/* when |x|<2**-14 */
    if(huge_f+x>one_f) return one_f+x;/* trigger inexact */
  }
  else k = 0;
  /* x is now in primary range */
  t  = x*x;
  if(k >= -125)
    SET_FLOAT_WORD(twopk,0x3f800000+(k<<23));
  else
    SET_FLOAT_WORD(twopk,0x3f800000+((k+100)<<23));
  c  = x - t*(P1_f+t*P2_f);
  if(k==0) 	return one_f-((x*c)/(c-(float)2.0f)-x);
  else 		y = one_f-((lo-(x*c)/((float)2.0f-c))-hi);
  if(k >= -125) {
    if(k==128) return y*2.0f*0x1p127f;
    return y*twopk;
  } else {
    return y*twopk*twom100_f;
  }
}
