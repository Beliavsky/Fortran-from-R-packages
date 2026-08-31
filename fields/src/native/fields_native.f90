module fields_native
use la_lapack, only: dpbtrf => pbtrf, dpbtrs => pbtrs
private
public :: bandsolve, bandchol, css, csstr, cvrf, dchold, ddfind, dlv, dmaket, drdfun, dsetup, &
          evlpoly, evlpoly2, expfn, ifind, igpoly, inpoly, inpoly2, mltdrb, mltdtd, &
          multwendlandg, wendlandfunction, multrb, radbas, radfun, rcss, rcssr, rcsswt, &
          rdist, rdist1, sortm
contains
!****fields is a package for analysis of spatial data written for
!****  # the R software environment .
!****  # Copyright (C) 2018
!****  # University Corporation for Atmospheric Research (UCAR)
!****  # Contact: Douglas Nychka, nychka@ucar.edu,
!****  # National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
!****  #
!****  # This program is free software; you can redistribute it and/or modify
!****  # it under the terms of the GNU General Public License as published by
!****  # the Free Software Foundation; either version 2 of the License, or
!****  # (at your option) any later version.
!****  # This program is distributed in the hope that it will be useful,
!****  # but WITHOUT ANY WARRANTY; without even the implied warranty of
!****  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!**** # GNU General Public License for more details.

!     wrapper for the lapack solve of a PD banded matrix using
!     the upper cholesky factor from bandchol
!
subroutine bandsolve( N, KD, NRHS, AB, LDAB, B, LDB, INFO )
integer N, KD, NRHS, LDAB, LDB, INFO
character UPLO
double precision AB(LDAB,*), B(LDB,*)
UPLO='U'
call DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
return
END

!     wrapper for the lapack cholesky of a PD banded matrix using
!     the upper cholesky factor from bandchol
!      
subroutine bandchol( N, KD, AB, LDAB, INFO)
integer N, KD, LDAB,INFO
character UPLO
double precision AB(LDAB,*)
UPLO='U'
call dpbtrf( UPLO, N, KD, AB, LDAB, INFO )
return
END


subroutine css(h,npoint,x,y,wght,sy,trace,diag,vlam, &
& ngrid,xg,yg,job,ideriv,ierr)

! COMPUTES A CUBIC SMOOTHING SPLINE FIT TO A SET OF DATA GIVEN  
! NPOINT(=NUMBER OF DATA VALUES) AND LAMBDA(=VALUE OF  
! THE SMOOTHING parameter). THE CODE IS ADAPTED FROM A  
! PROGRAM IN DEBOOR,C.(1978).A PRACTICAL GUIDE TO SPLINES.  
! SPRINGER-VERLAG : NEW YORK. AN O(NPOINT) ALGORITHM  
! SUGGESTED BY HUTCHINSON AND DEHOOG IS USED TO COMPUTE  
! LVRAGE VALUES AND CONSTRUCT CONFIDENCE INTERVALS.  
! Adapted from Randy Eubank Texas A&M , Arizona State 
!  
!  
! this subroutine solves the problem:  
!  
!   minimize  (1/n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + lambda*J(f)  
!    over f  
!   The solution will always be a piecewise cubic polynomial with join   
!   points at the x values. Natural boundary conditions are assumed: at the   
!   x(1) and x(npoints) the second and third derivatives of the slotuion   
!   will be zero   
!   All matrix calculations are done in double precision   
!  
!   Arguments of evss:   
!    h : natural log of lambda ( more convenient scale whepassing a   
!                               real*4)  
!        if h is passed with a value less than or equal -1000 no smoothing will   
!        be done and the spline will interploate the data points  
!    npoint: number of observations  
!    (x,y) : pairs of data points to be smoothed  
!    sy : on return predicted values of f at x  
!    wght : weights used in sum of squares. If the y have unequal   
!      variance then an appropriate choice for wght is the standard deviation  
!      (These weights are not normalized so multiplying by a constant   
!      will imply solving the minimization problem with a different smoothing   
!      parameter)    
!   trace: in matrix from Yhat= A(lambda)Y  with A(lambda) an nxn matrix  
!          trace= tr(A(lambda)) = " effective number of paramters"  
!   diag: diagonal elements of A(lambda) ( this is the mostt   
!         computationally intetnsive opertation in this subroutine.)  
!   vlam: value of the generalized cross-validation functyion (Used to  
!         select an appropriate value for lambda based on the data.)  
!   ngrid,xg,yg : on return, the ith deriv of the spline estimate is  
!                 evaluated on a grid, xg, with ngrid points.  The  
!                 values are returned in yg  
!  
!   ideriv: 0 = evaluate the function according to the job code  
!           1 = evaluate the first derivative of the function according  
!               to the job code  
!           2 = evaluate the second derivative of the function according  
!               to the job code  
!  
!   job: is a vector of three integers
!        (a,b,c)  (a=igcv, b=igrid, c=sorting)
!        a=0  evaluate spline at x values, return predicted values in sy  
!        a=1  same as a=0 plus returning values of trace, diag and vlam  
!        a=2  do no smoothing, interpolate the data  
!        a=3  same as a=1 but use the passed value invlam argument as  
!             a cost an offset factors in the diag vector  
!  
!        b=0  do not evaluate the spline at any grid points  
!        b=1  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
!             of the spline at the (unique) sorted,data points, xg, return yg  
!        b=2  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
!             of the spline at ngrid equally spaced points between x(1)  
!             and x(npoints)      
!        b=3  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
!             of the spline at ngrid points on grid passed in xg.
!             NOTE: extrapolation of the estimate beyond the range of   
!                     the x's will just be a linear function.  
!  
!
!
!         c=0    X's may not be in sorted order 
!         c=1    Assume that the X's are in sorted order
!         c=2 Do not sort the X's use the  current set of keys
!              Should only be used on a second call to smooth
!              same design
!
!  ierr : on return ierr>0 indicates all is not well :-(
!
parameter (NMAX=80000)
implicit double precision (a-h,o-z)
double precision h,trace,vlam
double precision wght(npoint),X(npoint),Y(npoint)
double precision sy(npoint),diag(npoint)
double precision xg(ngrid),yg(ngrid)
double precision A(NMAX,4),V(NMAX,7)
double precision P,SIXP,SIX1MP,cost
double precision ux(NMAX),uy(NMAX), uw(NMAX),ud(NMAX),utr
integer imx(NMAX)
integer idx(NMAX)
integer npoint,igcv,igrid, isort,  job(3)

eps= 1d-7
cost=1.0
offset= 0
igcv= job(1)
igrid=job(2)
isort=job(3)

!  
if( npoint.gt.NMAX) then
ierr=1
return
endif
! initialize unique vector and two pointers

do 5 k=1,npoint
ux(k)=x(k)

idx(k)=k
5 continue

! sort vector X along with the keys in imx
!
!       initialize the indices  imx
!
if( isort.le.1) then

do 6 k=1,npoint
imx(k)=k
6 continue
endif
!   
!    sort the X's permuting the indices
! 
if(isort.eq.0) then
call sortm( ux, imx,npoint)
!   the rank of the value x( imx(k)) is k
endif
! put y and the weights in the  sorted order 
!  
!**** construct vector consisting of unique observations.  

jj=1
ind= imx(1)
ux(jj)= x(ind)
uy(jj)= y(ind)
uw(jj)= wght(ind)
idx(1)=jj
!      normalize eps by the range of the X values
eps= eps/(  x( imx(npoint)) - x( imx(1)) )
!
!
do 10 k=2,npoint
!   we are looping through ranks but this is not how the 
!   order of the X are stored. The location of the kth smallest
!   is at idx(k) 
kshuf= imx(k)
if( abs( ux(jj)-x(kshuf)).lt.eps) then
!**** we have a repeat observation, update weight and the weighted   
!***** average at this point  
temp1=  1.0d0/uw(jj)**2
temp2=  1.0d0/wght(kshuf)**2
temp3 = (temp1 + temp2)
uy(jj)=  ( uy(jj)*temp1 + y(kshuf)*temp2)/temp3
uw(jj)= 1.0d0/dsqrt(temp3)
else
jj= jj+1
ux(jj)= x(kshuf)
uy(jj)= y(kshuf)
uw(jj)= wght(kshuf)
endif
!  save the value that indexes unique values to repliacted ones.
!    x(k) corresponds to the unique X at idx(k)  
idx(k)=jj
10 continue
nunq= jj

itp=0
if(igcv.eq.2) itp=1

!   handle condition for interpolation        

if(itp.eq.0) then
P=1.d0/(npoint*dexp(h)+ 1.d0)
else
P=1.d0
endif


call dsetup(uX,uW,uY,nunq,V,A(1,4),NMAX,itp,ierr)
!****  check for duplicate X's if so exit   
if(ierr.gt.0) then
return
endif

call dchold(P,V,A(1,4),nunq,A(1,3),A(1,1),NMAX)

! compute predicted values   

SIX1MP=6.d0*(1.d0-P)
if(itp.eq.0) then
DO 61 I=1,Nunq
a(i,1)=uY(I) - SIX1MP*(uW(I)**2)*A(I,1)
61 continue

! fill in predicted values taking into account repeated data  
do 70 k=1,npoint
jj= idx(k)
sytemp= a(jj,1)
! 
! unscramble the smoothed Y's
kshuf= imx(k)

sy(kshuf)=sytemp
70 continue
else
do 60 i=1,nunq
a(i,1)=uy(i)
60 continue
endif

! return estimates on unique x's if igrid =1  
if(igrid.eq.1) then
do 65 i=1,nunq
xg(i)=ux(i)
yg(i)=a(i,1)
65 continue
ngrid=nunq
endif
if(igrid.ge.1) then
!  
!********* evaluate spline on grid  
! piecewise cubic COEFFICIENTS ARE STORED IN A(.,2-4).   

SIXP=6.d0*P
DO 62 I=1,nunq
A(I,3)=A(I,3)*SIXP
62 continue
NPM1=nunq - 1
DO 63 I=1,NPM1
A(I,4)=(A(I+1,3)-A(I,3))/V(I,4)
A(I,2)=(A(I+1,1)-A(I,1))/V(I,4) &
& -(A(I,3)+A(I,4)/3.*V(I,4))/2.*V(I,4)
63 continue
!  
!  create equally spaced x's if asked for ( igrid=2)  
!  
if( igrid.eq.2) then
step= (ux(nunq)-ux(1))/(ngrid-1)
xg(1)=ux(1)
do 190 j=2,ngrid-1
xg(j)= xg(j-1)+step
190 continue
xg(ngrid)=ux(nunq)
endif


uxmin= ux(1)
uxmax= ux(nunq)

a1= a(1,1)
an= a(nunq,1)

b1= a(1,2)

d= ux(nunq)- ux(nunq-1)
ind= nunq-1
bn= a(ind,2) + a(ind,3)*d + a(ind,4)*(d**2)/2.d0

! evalute spline by finding the interval containing the evaluation   
! point and applying the cubic formula for the curve estiamte  
! finding the interval such that ux(ind) <=xg(j) < ux( ind+1)  
! is done using a bisection search   

do 195 j=1,ngrid
lint= ifind(xg(j),ux,nunq)
if( lint.eq.0) then
d= xg(j)-uxmin

if (ideriv .eq. 0) &
& yg(j)= a1 + b1*d
if (ideriv .eq. 1) &
& yg(j)=  b1
if (ideriv .eq. 2) &
& yg(j)= 0.0
endif
if( lint.eq.nunq) then
d= xg(j)-uxmax

if (ideriv .eq. 0) &
& yg(j)= an + bn*d
if (ideriv .eq. 1) &
& yg(j)=  bn
if (ideriv .eq. 2) &
& yg(j)= 0.0
endif
if( ((lint.ge.1 ) .and.( lint.lt.nunq))) then
ind=lint
!            a1=a(ind,1)  
!            a2=a(ind,2)  
!            b=a(ind,3)/2.d0  
!            c=a(ind,4)/6.d0  
!  
d= xg(j)-ux(ind)

if (ideriv .eq. 0) &
& yg(j)= a(ind,1) + a(ind,2)*d + a(ind,3)*(d**2)/2.d0 &
& + a(ind,4)*(d**3)/6.d0
if (ideriv .eq. 1) &
& yg(j)= a(ind,2) + a(ind,3)*d + a(ind,4)*(d**2)/2.d0
if (ideriv .eq. 2) &
& yg(j)= a(ind,3) + a(ind,4)*d
endif
!  
195 continue
endif
!****end of evaluation block  

if((igcv.eq.1).or.( igcv.eq.3)) then
if( igcv.eq.3)  then
cost= diag(1)
offset=diag(2)
endif
!*****                       begin computing gcv and trace  
! COMPUTE LVRAGE VALUES ,THE VARIANCE ESTIMATE   
!     SGHAT2 AND CONFIDENCE INTERVALS.  
!  
call dLV(nunq,V,uw,SIX1MP,utr,ud,NMAX)

rss=0.d0
trace=0.d0
vlam=0.d0

do 100 i=1,nunq
!            rss= rss + ((uy(i)-a(i,1))/uw(i))**2  
trace= trace +ud(i)
100 continue

do 110 k=1,npoint
kshuf= imx(k)
jj= idx(k)
diag(kshuf)= ud(jj)
rss= rss + ( (y(kshuf)- sy(kshuf))/wght(kshuf) )**2
110 continue
ctrace=  2+  cost*( trace-2)
if( (npoint -ctrace -offset) .gt. 0.d0) then
vlam= (rss/npoint)/( 1- (ctrace-offset)/npoint)**2
else
vlam=1e20
endif

endif



return

END
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html


subroutine csstr(h,nobs,x,y,wght,c,offset,trace,vlam,work,ierr)
parameter(mxM=20000)
implicit double precision (a-h,o-z)
double precision h,trace, vlam,c,offset
double precision x(nobs),y(nobs),wght(nobs)
double precision work(nobs),diag(mxM),dumm1(1),dumm2(1)
integer job(3),ideriv,ierr, ndum
data ideriv/0/
job(1)=3
job(2)=0
job(3)=0
diag(1)=c
diag(2)=offset
ndum=1
call css(h,nobs,x,y,wght,work,trace,diag,vlam,ndum,dumm1,dumm2, &
& job,ideriv,ierr)

return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html



!** finds value of h minimizing the  generalized cross-validation

subroutine cvrf &
& (h,nobs,x,y,wt,sy,trace,diag,din,dout,cv,ierr)
implicit double precision (a-h,o-z)
double precision h,trace,cv
double precision x(nobs),y(nobs),wt(nobs)
double precision sy(nobs),diag(nobs),dumm1(1),dumm2(1)
double precision din(10), dout(10)
integer ngrid, ierr, job(3),ideriv
data job/3,0,0/
data ideriv,ngrid/0,1/
call rcss(h,nobs,x,y,wt,sy,trace,diag,cv, &
& ngrid,dumm1,dumm2,job,ideriv,din,dout,ierr)
nit= int( dout(1))
trace=dout(3)

!
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html





SUBROUTINE dchold(P,V,QTY,NPOINT,U,QU,NMAX)
! CONSTRUCTS THE UPPER THREE DIAGONALS IN V(I,J),I=2,
!   NPOINT-1,J=1,3 OF THE MATRIX 6*(1-P)*Q-TRANSP*
!   (D**2)*Q + P*R . (HERE R IS THE MATRIX PROPORTIONAL
!   TO Q-TRANSP*KSUBN*Q , WHERE KSUBN IS THE MATRIX 
!   WITH ELEMENTS K(X(I),X(J)) AND K IS THE USUAL
!   KERNEL ASSOCIATED WITH CUBIC SPLINE SMOOTHING.)
!   THE CHOLESKY DECOMPOSITION OF THIS MATRIX IS COMPUTED
!   AND STORED IN V(.,1-3) AND THE EQUATION SYSTEM FOR 
!   THE QUADRATIC COEFFICIENTS OF THE SPLINE IN ITS 
!   PIECEWISE POLYNOMIAL REPRESENTATION IS SOLVED . THE 
!   SOLUTION IS RETURNED IN U.
!
double precision P,QTY(NMAX),QU(NMAX),U(NMAX),V(NMAX,7)
double precision SIX1MP,TWOP,RATIO,PREV
INTEGER NPOINT,I,NPM1,NPM2
!
NPM1=NPOINT - 1
!****  CONSTRUCT 6*(1-P)*Q-TRANSP.*(D**2)*Q + P*R
SIX1MP=6.d0*(1.d0 - P)
TWOP=2.d0*P
DO 2 I=2,NPM1
V(I,1)=SIX1MP*V(I,5) + TWOP*(V(I-1,4)+V(I,4))
V(I,2)=SIX1MP*V(I,6) + P*V(I,4)
V(I,3)=SIX1MP*V(I,7)
2 continue
NPM2=NPOINT - 2
IF(NPM2 .GE. 2)GO TO 10
U(1)=0.d0
U(2)=QTY(2)/V(2,1)
U(3)=0.d0
GO TO 41
!  FACTORIZATION : FACTORIZE THE MATRIX AS L*B-INV*
!     L-TRANSP , WHERE B IS A DIAGONAL MATRIX AND L
!     IS UPPER TRIANGULAR.
10 DO 20 I=2,NPM2
RATIO=V(I,2)/V(I,1)
V(I+1,1)=V(I+1,1) - RATIO*V(I,2)
V(I+1,2)=V(I+1,2) - RATIO*V(I,3)
V(I,2)=RATIO
RATIO=V(I,3)/V(I,1)
V(I+2,1)=V(I+2,1) - RATIO*V(I,3)
V(I,3)=RATIO
20 continue
!  FORWARD SUBSTITUTION
U(1)=0.d0
V(1,3)=0.d0
U(2)=QTY(2)
DO 30 I=2,NPM2
U(I+1)=QTY(I+1) - V(I,2)*U(I) -V(I-1,3)*U(I-1)
30 continue
!  BACK SUBSTITUTION
U(NPOINT)=0.d0
U(NPM1)=U(NPM1)/V(NPM1,1)
DO 40 I=NPM2,2,-1
U(I)=U(I)/V(I,1) - U(I+1)*V(I,2) - U(I+2)*V(I,3)
40 continue
!  CONSTRUCT Q*U
41 PREV=0.d0
DO 50 I=2,NPOINT
QU(I)=(U(I)-U(I-1))/V(I-1,4)
QU(I-1)=QU(I) - PREV
PREV=QU(I)
50 continue
QU(NPOINT)=-QU(NPOINT)
!
RETURN
END
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

!
subroutine ddfind( nd,x1,n1, x2,n2, D0,ind,rd,Nmax, iflag)

integer nd,n1,n2, ind(nmax,2)
integer kk, i,j, ic

double precision x1(n1,nd), x2(n2,nd), D0, rd(Nmax), D02, dtemp

!****   counter  for accumulating close points
kk=0
D02= D0**2
do  15 i= 1, n1
do 10 j =1,n2

!
!** accumulate squared differences
!
dtemp= 0.0
do 5 ic= 1, nd
dtemp= dtemp + (x1(i,ic) - x2(j,ic))**2
if( dtemp.gt.D02) goto 10
5 continue

!****       dtemp is less than D0 so save it as a close point

kk=kk+1

!**** check if there is still space 
if( kk .gt. Nmax) then
iflag= -1
goto 20
else
ind(kk,1)= i
ind(kk,2)= j
rd(kk)= sqrt( dtemp)
endif


10 continue
15 continue

Nmax=kk
20 continue


return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html





SUBROUTINE dLV(NPOINT,V,WGHT,SIX1MP,TR,LEV,NMAX)
! CONSTRUCTS THE UPPER THREE DIAGONALS OF (6*(1-P)*
!   Q-TRANSP*(D**2)*Q + P*R)-INV USING THE RECURSION
!   FORMULA IN HUTCHINSON,M.F. AND DEHOOG,F.R.(1985).
!   NUMER. MATH. 47,99-106, AND STORES THEM IN V(.,5-7).
!   THESE ARE USED IN POST AND PRE MULTIPLICATION BY
!   Q-TRANSP AND Q TO OBTAIN THE DIAGONAL ELEMENTS OF
!   THE HAT MATRIX WHICH ARE STORED IN THE VECTOR LEV.
!     THE TRACE OF THE HAT MATRIX IS RETURNED IN TR.
!
double precision V(NMAX,7),TR,W1,W2,W3,SIX1MP
double precision wght(NMAX)
double precision LEV(npoint)
INTEGER NPM1,NPM2,NPM3,NPOINT
!
NPM1=NPOINT - 1
NPM2=NPOINT - 2
NPM3=NPOINT - 3
!   RECURSION FOR DIAGONALS OF INVERSE MATRIX
V(NPM1,5)=1/V(NPM1,1)
V(NPM2,6)=-V(NPM2,2)*V(NPM1,5)
V(NPM2,5)=(1/V(NPM2,1)) - V(NPM2,6)*V(NPM2,2)
DO 10 I=NPM3,2,-1
V(I,7)=-V(I,2)*V(I+1,6) - V(I,3)*V(I+2,5)
V(I,6)=-V(I,2)*V(I+1,5) - V(I,3)*V(I+1,6)
V(I,5)=(1/V(I,1))- V(I,2)*V(I,6) - V(I,3)*V(I,7)
10 CONTINUE
!   POSTMULTIPLY BY (D**2)*Q-TRANSP AND PREMULTIPLY BY Q TO
!     OBTAIN DIAGONALS OF MATRIX PROPORTIONAL TO THE 
!     IDENTITY MINUS THE HAT MATRIX.
W1=1.d0/V(1,4)
W2= -1.d0/V(2,4) - 1.d0/V(1,4)
W3=1.d0/V(2,4)
V(1,1)=V(2,5)*W1
V(2,1)=W2*V(2,5) + W3*V(2,6)
V(2,2)=W2*V(2,6) + W3*V(3,5)
LEV(1)=1.d0 - (WGHT(1)**2)*SIX1MP*W1*V(1,1)
LEV(2)=1.d0 - (WGHT(2)**2)*SIX1MP*(W2*V(2,1) + W3*V(2,2))
TR=LEV(1) + LEV(2)
DO 20 I=4,NPM1
W1=1.d0/V(I-2,4)
W2= -1.d0/V(I-1,4) - 1.d0/V(I-2,4)
W3=1.d0/V(I-1,4)
V(I-1,1)=V(I-2,5)*W1 + V(I-2,6)*W2 + V(I-2,7)*W3
V(I-1,2)=V(I-2,6)*W1 + V(I-1,5)*W2 + V(I-1,6)*W3
V(I-1,3)=V(I-2,7)*W1 + V(I-1,6)*W2 + V(I,5)*W3
LEV(I-1)=1.d0 - (WGHT(I-1)**2)*SIX1MP*(W1*V(I-1,1) &
& + W2*V(I-1,2) + W3*V(I-1,3))
TR= TR + LEV(I-1)
20 CONTINUE
W1=1.d0/V(NPM2,4)
W2= -1.d0/V(NPM1,4) - 1.d0/V(NPM2,4)
W3=1.d0/V(NPM1,4)
V(NPOINT,1)=V(NPM1,5)*W3
V(NPM1,1)=V(NPM2,5)*W1 + V(NPM2,6)*W2
V(NPM1,2)=V(NPM2,6)*W1 + V(NPM1,5)*W2
LEV(NPM1)=1.d0 - (WGHT(NPM1)**2)*SIX1MP*(W1*V(NPM1,1) &
& + W2*V(NPM1,2))
LEV(NPOINT)=1.d0 - (WGHT(NPOINT)**2)*SIX1MP*W3*V(NPOINT,1)
TR= TR + LEV(NPM1) + LEV(NPOINT)
RETURN
END
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html


subroutine dmaket(m,n,dim,des,lddes,npoly,t,ldt, &
& wptr,info,ptab,ldptab)
integer m,n,dim,lddes,npoly,ldt,wptr(dim),info,ptab(ldptab,dim)
double precision des(lddes,dim),t(ldt,npoly)
!
! Purpose: create t matrix and append s1 to it.
!
! On Entry:
!   m			order of the derivatives in the penalty
!   n			number of rows in des
!   dim			number of columns in des
!   des(lddes,dim)	variables to be splined
!   lddes		leading dimension of des as declared in the
!			calling program
!   ldt			leading dimension of t as declared in the
!			calling program
!
!   npoly		dimension of polynomial part of spline
! On Exit:
!   t(ldt,npoly+ncov1)	[t:s1]
!   info 		error indication
!   			   0 : successful completion
!		 	   1 : error in creation of t
! Work Arrays:
!   wptr(dim)		integer work vector
!
! Subprograms Called Directly:
!	Other - mkpoly
!
!
integer i,j,k,kk,tt,nt,bptr,eptr
!
info = 0
!      npoly = mkpoly(m,dim)
do 5 j=1,n
t(j,1)=1.0
5 continue
nt = 1
if (npoly .gt. 1) then
do 10 j=1,dim
nt = j + 1
wptr(j) = nt
ptab(nt,j)= ptab(nt,j) +1
do 15 kk = 1, n
t(kk,nt)= des(kk,j)
15 continue
!             call dcopy(n,des(1,j),1,t(1,nt),1)
10 continue
!
!     get cross products of x's in null space for m>2
!
!     WARNING: do NOT change next do loop unless you fully understand:
!              This first gets x1*x1, x1*x2, x1*x3, then
!              x2*x2, x2*x3, and finally x3*x3 for dim=3,n=3
!              wptr(1) is always at the beginning of the current
!	       level of cross products, hence the end of the
!	       previous level which is used for the next.
!	       wptr(j) is at the start of xj * (previous level)
!
do 50 k=2,m-1
do 40 j=1,dim
bptr = wptr(j)
wptr(j) = nt + 1
eptr = wptr(1) - 1
do 30 tt=bptr,eptr
nt = nt + 1
do 21 jj= 1,dim
ptab(nt,jj)= ptab(tt,jj)
21 continue
ptab( nt,j)= 1+ ptab( nt,j)
do 20 i=1,n
t(i,nt) = des(i,j) * t(i,tt)
20 continue
30 continue
40 continue
50 continue
if (nt .ne. npoly) then
info = 1
return
endif
endif
!
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

subroutine drdfun(n, d2, par)
double precision d2(n), par(2), dtemp
integer n
if( int(par(2)).eq.0) then

do 5 k =1,n
d2(k)= par(1)*(d2(k))**( par(1)-1)
5 continue
else
do 6 k=1,n
dtemp= d2(k)
if( dtemp.GE.1e-35)  then
!
! NOTE factor of 2 adjusts for log being applied to 
! distance rather than squared distance
d2(k)=  (par(1)*log(dtemp) +1)*(dtemp)**( par(1)-1)/2
else
d2(k)=0.0
endif
6 continue
endif
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

SUBROUTINE dsetup(X,WGHT,Y,NPOINT,V,QTY,NMAX,itp,ierr)
!   PUT DELX=X(.+1)-X(.) INTO V(.,4)
!   PUT THE THREE BANDS OF THE MATRIX Q-TRANSP*D INTO
!     V(.,1-3)
!   PUT THE THREE BANDS OF (D*Q)-TRANSP*(D*Q) AT AND 
!     ABOVE THE DIAGONAL INTO V(.,5-7)
!   HERE Q IS THE TRIDIAGONAL MATRIX OF ORDER (NPOINT
!     -2,NPOINT) THAT SATISFIES Q-TRANSP*T=0 AND WGHT
!     IS THE DIAGONAL MATRIX WHOSE DIAGONAL ENTRIES 
!     ARE THE SQUARE ROOTS OF THE WEIGHTS USED IN THE 
!     PENALIZED LEAST-SQUARES CRITERION
!
implicit double precision (a-h,o-z)
double precision WGHT(NMAX),X(NMAX),y(NMAX)
double precision QTY(NMAX),V(NMAX,7)
double precision DIFF,PREV
INTEGER NPOINT,I,NPM1
!
NPM1=NPOINT -1
V(1,4)=X(2)-X(1)
if(v(1,4).eq.0.d0) then
ierr=5
return
endif

DO 11 I=2,NPM1
V(I,4)=X(I+1) - X(I)
if(v(I,4).eq.0.d0) then
ierr=5
return
endif
if(itp.eq.0) then
V(I,1)=WGHT(I-1)/V(I-1,4)
V(I,2)=-WGHT(I)/V(I,4) - WGHT(I)/V(I-1,4)
V(I,3)=WGHT(I+1)/V(I,4)
else
V(I,1)=1.d0/V(I-1,4)
V(I,2)=-1.d0/V(I,4) - 1.0/V(I-1,4)
V(I,3)=1.d0/V(I,4)
endif
11 continue
!
V(NPOINT,1)=0.d0
DO 12 I=2,NPM1
V(I,5)=V(I,1)**2 + V(I,2)**2 + V(I,3)**2
12 continue
IF(NPM1 .LT. 3)GO TO 14
DO 13 I=3,NPM1
V(I-1,6)=V(I-1,2)*V(I,1) + V(I-1,3)*V(I,2)
13 continue
14 V(NPM1,6)=0.d0
IF(NPM1 .LT. 4)GO TO 16
DO 15 I=4,NPM1
V(I-2,7)=V(I-2,3)*V(I,1)
15 continue
16 V(NPM1-1,7)=0.d0
V(NPM1,7)=0.d0
!
!  CONSTRUCT Q-TRANSP. * Y IN QTY
PREV=(Y(2) - Y(1))/V(1,4)
DO 21 I=2,NPM1
DIFF=(Y(I+1)-Y(I))/V(I,4)
QTY(I)=DIFF - PREV
PREV=DIFF
21 continue
!
RETURN
END
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html



subroutine evlpoly(x,n,coef,j,result)

! evaluates a polynomial: coef(1) + sum_i= 2,j coef(i)x**i

integer j,n, i
double precision x(n), result(n), coef(j)
double precision temp, tempx, temp2

do 10 i = 1,n

temp= coef(1)
tempx= x(i)
temp2= tempx

!      temp is set to constant now loop over powers

do 20 kk= 2, j
temp= temp + coef(kk)* temp2
temp2= temp2*tempx
20 continue

result(i)= temp

10 continue

return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html



subroutine evlpoly2(x,n,nd,ptab,j,coef,result)

! evaluates a polynomial: coef(1) + sum_i= 2,j coef(i)x**i

integer j,n, i, nd
double precision x(n,nd), result(n), coef(j)
integer ptab(j,nd)
double precision temp,  temp2

do 10 i = 1,n

temp= 0

! for a given vector accumlate the polynomial terms

do 20 kk= 1, j
temp2 =1.0

do 30 l=1, nd
if( ptab(kk,l).ne.0) then
temp2= temp2* (x(i,l)**ptab(kk,l))
endif
30 continue

temp = temp + temp2*coef(kk)

20 continue

result(i)= temp

10 continue

return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html


subroutine expfn(n,d2, par)
double precision d2(n), par(1)
integer n

do 5 k =1,n
d2(k)=  exp(-1*d2(k)**(par( 1)/2))
5 continue

return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html


INTEGER FUNCTION IFIND(X,XK,N)
!  FIND I SUCH THAT XK(I) LE X LT XK(I+1)                              
!  IFIND=0  IF X LT XK(1)                                              
!  IFIND=N  IF X GT XK(N)                                              
!  J F MONAHAN  JAN 1982  DEPT OF STAT, N C S U, RALEIGH, NC 27650     
double precision X,XK(N)
IFIND=0
IF(X.LT.XK(1)) RETURN
IFIND=N
IF(X.GE.XK(N)) RETURN
IL=1
IU=N
1 IF(IU-IL.LE.1) GO TO 4
I=(IU+IL)/2
!      IF(X-XK(I)) 2,5,3                                                
IF( (X-XK(I)).eq.0) go to 5
IF( (X-XK(I)).gt.0)  go to 3
! following used to have line number 2      
IU=I
GO TO 1
3 IL=I
GO TO 1
4 IFIND=IL
RETURN
5 IFIND=I
RETURN
END
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html


subroutine igpoly(nx, xg,ny,yg,np,xp,yp,ind)
!----------------------------------------------------------------------
!     This subroutine determines whether or not an 2-d point (xg(i),yg(j))
!     element is inside a (closed) set of points xp,yp that
!     are assumed to form polygon.
!     result is an  matrix indicating comparision for all 
!     combinations of xg and yg
!----------------------------------------------------------------------

integer np
!    # of points in polygon
integer nx,ny
!     # points to check 
real xg(nx)
!     x grid values to check
real yg(ny)
!     y grid values to check 
real    xp(np)
!     2d-locations of polygon
real  yp(np)
real  x1, x2, y1,y2
!     min and max of x and y
real  temp, xt, yt
integer ind(nx,ny)
! THE ANSWER : ind(i)=1 if point xd(i),yd(i) is 
!     in polygon 0 otherwise 
integer in
x1= xp(1)
x2= xp(2)
y1= yp(1)
y2= yp(1)
!
!     find the minima and maxima of the polygon coordinates
!     i.e. the smallest rectangle containing the polygon.
do j = 1,np

temp= xp(j)
if( temp.lt.x1)  then
x1 = temp
endif
if( temp.gt.x2) then
x2 = temp
endif

temp= yp(j)
if( temp.lt.y1)  then
y1 = temp
endif
if( temp.gt.y2) then
y2 = temp
endif
enddo

! loop over all 
! grid points

do i = 1, nx
do j = 1,ny
xt= xg(i)
yt= yg(j)
! quick test that point is inside the bounding rectangle
! of the polygon

if( (xt.le. x2).and. (xt.ge.x1).and. &
& (yt.le.y2).and.(yt.ge.y1) ) then
call inpoly2(xt,yt,np,xp,yp, in)
ind(i,j)=in
else
ind(i,j)= 0
endif

enddo
enddo

return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html


subroutine inpoly(nd, xd,yd,np,xp,yp,ind)
!----------------------------------------------------------------------
!     This subroutine determines whether or not an 2-d point (xd(j),yd(j))
!     element is inside a (closed) set of points xp,yp that
!     are assumed to form polygon.
!----------------------------------------------------------------------

integer np
!     # of points in polygon
integer nd
!     # points to check 
real xd(nd)
!     2d-locations to check
real yd(nd)
real    xp(np)
!     2d-locations of polygon
real  yp(np)
real  x1, x2, y1,y2
!     min and max of x and y
real  temp, xt, yt
integer ind(nd)
!     THE ANSWER : ind(i)=1 if point xd(i),yd(i) is 
!  in polygon 0 otherwise 
integer in
x1= xp(1)
x2= xp(2)
y1= yp(1)
y2= yp(1)
!
!     find the minima and maxima of the polygon coordinates
!     i.e. the smallest rectangle containing the polygon.
do j = 1,np

temp= xp(j)
if( temp.lt.x1)  then
x1 = temp
endif
if( temp.gt.x2) then
x2 = temp
endif

temp= yp(j)
if( temp.lt.y1)  then
y1 = temp
endif
if( temp.gt.y2) then
y2 = temp
endif
enddo

do j = 1,nd
xt= xd(j)
yt= yd(j)
! quick test that point is inside the bounding rectangle
! if not it is not inside polygon
if( (xt.le. x2).and. (xt.ge.x1).and. &
& (yt.le.y2).and.(yt.ge.y1) ) then
call inpoly2(xt,yt,np,xp,yp, in)
ind(j)=in
else
ind(j)= 0
endif

enddo

return
end

subroutine inpoly2(xpnt,ypnt,np,xp,yp,in)
!
parameter (pi=3.14159265358979,ttpi=2.*pi)
!
dimension xp(np),yp(np)
real xpnt, ypnt
integer in
!
!----------------------------------------------------------------------
!
!  THE VALUE OF THIS FUNCTION IS NON-ZERO IF AND ONLY IF (XPNT,YPNT) 
!  IS INSIDE OR *ON* THE POLYGON DEFINED BY THE POINTS (XP(I),YP(I)), 
!  FOR I FROM 1 TO NP.
!
!  THE INPUT POLYGON DOES NOT HAVE TO BE A CLOSED POLYGON, THAT IS
!  IT DOES NOT HAVE TO BE THAT (XP(1),YP(1) = (XP(NP,YP(NP)).
!
!----------------------------------------------------------------------
!
!  DETERMINE THE NUMBER OF POINTS TO LOOK AT (DEPENDING ON WHETHER THE
!  CALLER MADE THE LAST POINT A DUPLICATE OF THE FIRST OR NOT).
!
if (xp(np).eq.xp(1) .and. yp(np).eq.yp(1)) then
npts = np-1
else
npts = np
end if

in = 0
!     ASSUME POINT IS OUTSIDE

! --- ------------------------------------------------------------------
! --- CHECK TO SEE IF THE POINT IS ON THE POLYGON.
! --- ------------------------------------------------------------------

do ipnt = 1,npts
if (xpnt .eq. xp(ipnt) .and. ypnt .eq. yp(ipnt) ) then
in = 1
goto 999
! EARLY EXIT
endif
enddo

! --- ------------------------------------------------------------------
! --- COMPUTE THE TOTAL ANGULAR CHANGE DESCRIBED BY A RAY EMANATING 
! --- FROM THE POINT (XPNT,YPNT) AND PASSING THROUGH A POINT THAT 
! --- MOVES AROUND THE POLYGON.
! --- ------------------------------------------------------------------

anch = 0.

inxt = npts
xnxt = xp(npts)
ynxt = yp(npts)
anxt = atan2(ynxt-ypnt, xnxt-xpnt)

do 100 ipnt=1,npts

ilst = inxt
xlst = xnxt
ylst = ynxt
alst = anxt

inxt = ipnt
xnxt = xp(inxt)
ynxt = yp(inxt)
anxt = atan2(ynxt-ypnt, xnxt-xpnt)

adif = anxt-alst

if (abs(adif) .gt. pi) adif = adif - sign(ttpi,adif)

anch = anch + adif

100 continue

! --- ------------------------------------------------------------------
! --- IF THE POINT IS OUTSIDE THE POLYGON, THE TOTAL ANGULAR CHANGE 
! --- SHOULD BE EXACTLY ZERO, WHILE IF THE POINT IS INSIDE THE POLYGON,
! --- THE TOTAL ANGULAR CHANGE SHOULD BE EXACTLY PLUS OR MINUS TWO PI.
! --- WE JUST TEST FOR THE ABSOLUTE VALUE OF THE CHANGE BEING LESS 
! --- THAN OR EQUAL TO PI.
! --- ------------------------------------------------------------------

if (abs(anch) .ge. pi) in = 1

999 continue
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html



!** evaluates partial derivatives of radial basis functions with
!** nodes at x2 and at the points x1
!
subroutine mltdrb( nd,x1,n1, x2, n2, par, c,h,work)
implicit double precision (a-h,o-z)
integer nd,n1,n2,ic, ivar, ir, j
double precision par(2), x1(n1,nd)
double precision  x2(n2,nd), c(n2), h(n1, nd)
double precision sum, sum2
double precision work( n2)
!      double precision ddot
do 1000 ivar=1, nd
!****** work aray must be dimensioned to size n2
! **** loop through columns of output matrix K
!***  outermost loop over columns of x1 and x2 should
!***  help to access adjacent values in memory.          
do 5 ir= 1, n1
! evaluate all basis functions at  x1(j,.)       
do 10 j =1,n2
!  zero out sum accumulator
sum=0.0
do 15  ic=1,nd
!** accumulate squared differences
sum= sum+ (x1(ir,ic)- x2(j,ic))**2
15 continue
work(j)=sum
10 continue
!**** evaluate squared distances  with basis functions. 
call drdfun( n2,work(1) ,par(1) )
do 11 j= 1, n2
work( j)= 2.0*work(j)*(x1(ir,ivar)- x2(j,ivar))
11 continue
!*****now the dot product you have all been waiting for!
sum2= 0.0
do 12 j = 1, n2
sum2 = sum2  + work(j)*c(j)
12 continue
!     h(ir,ivar)= ddot( n2, work(1), 1, c(1),1)
h(ir,ivar) = sum2
5 continue
1000 continue
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html



subroutine mltdtd( nd,x1,n1,np,ptab, d,h)
implicit double precision (a-h,o-z)
integer nd,n1,np, ivar

double precision x1(n1,nd)
double precision   d(np), h(n1, nd)
double precision work,  prod, xp,xx
integer ptab(np,nd)
!  outer most loop is over the variables w/r partial derivative

do 1000 ivar=1, nd
!****** work aray must be dimensioned to size np

do 5 ir= 1, n1
! next loop is over rows of x1

!

! evaluate all partials of polynomials  at  x1(j,.)       
! take ddot product of this vector with d 
! this is the element to return in h(ir,ivar)
work=0.0
do 10 j =1,np
prod=0.0
ipv= ptab( j,ivar)
if( ipv.gt.0) then
prod=1.0
do 11 k= 1, nd
ip= ptab(j,k)
! ip is the power of the kth variable in the jth poly
if( ip.eq.0) goto 11
xx= x1(ir,k)
!**
if( k.eq.ivar) then
if( ip.eq.1) then
xp=1.0
else
xp= (ip)* xx**(ip-1)
endif
else
xp= xx**(ip)
endif
prod=prod*xp
11 continue
endif
work= work + prod* d(j)
10 continue

!
h(ir,ivar)=work
5 continue
1000 continue
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

subroutine multWendlandG( mx, my, deltaX, deltaY, &
& nc, center, coef, h, flag)
integer mx, my, nc, flag
double precision  deltaX, deltaY, center(nc,2), coef(nc)
double precision h(mx,my)
! 
integer j, k, l, m1, m2, n1, n2
double precision  kstar, lstar, d
do j = 1, nc
kStar= center(j,1)
lStar= center(j,2)
m1 =  max( ceiling(-deltaX + kStar),  1)
m2 =  min(   floor( deltaX + kStar), mx)
n1 =  max( ceiling(-deltaY + lStar),  1)
n2 =  min(   floor( deltaY + lStar), my)
do l = n1, n2
do k = m1, m2
d = dsqrt( ((k-kStar)/deltaX)**2 + ((l- lStar)/deltaY)**2)
h(k,l) = h(k,l) + wendlandFunction( d)* coef(j)
!             h(k,l) = h(k,l) + 2
enddo
enddo
enddo
flag=0
return
end

double precision function wendlandFunction(d)
double precision d
if( d.GE.1) then
wendlandFunction = 0
else
wendlandFunction = ((1-d)**6) * (35*d**2 + 18*d + 3)/3
endif
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

!** evaluates radial basis functions 
!**** K_ij= radfun( distance( x1_i, x2_j))
!**** and does the multplication  h= Kc
!****  K is n1Xn2
!****  h is n1Xn3
!***** c is n2xn3


subroutine multrb( nd,x1,n1, x2,n2, par, c,n3,h,work)
implicit double precision (a-h,o-z)
integer nd,n1,n2,n3,ic, jc,j
double precision par(2),x1(n1,nd), x2(n2,nd), c(n2,n3), h(n1,n3)
double precision work( n2)
!       double precision ddot
double precision sum, sum2
! radfun has an explicit module interface.

!****** work aray must be dimensioned to size n1
! **** loop through columns of output matrix K
do 5 ir= 1, n1
do 10 j =1,n2
sum=0.0
do 15  ic=1,nd
!** accumulate squared differences
sum= sum+ (x1(ir,ic)- x2(j,ic))**2
15 continue
work(j)=radfun( sum, par(1), par(2))
10 continue
!***** dot product for matrix multiplication
do 30 jc=1,n3
sum2= 0.0
do 12 j = 1, n2
sum2 = sum2  + work(j)*c(j, jc)
12 continue
h(ir,jc) = sum2
!      h(ir,jc)= ddot( n2, work, 1, c(1,jc),1)
30 continue
5 continue
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

!**** subroutine to fill in the omega ( or K) matrix for 
!**** ridge regression S funcion
!**** K_ij= radfun( distance( x1_i, x2_j))
!
subroutine radbas( nd,x1,n1, x2,n2, par, k)
integer nd,n1,n2,ic
double precision par(2), x1(n1,nd), x2(n2,nd), k(n1,n2)
double precision xtemp
! **** loop through columns of output matrix K
!*** outer most loop over columns of x1 and x2 should reduce memory swaps 
do  ic= 1, nd
do  j =1,n2
xtemp= x2(j,ic)
do   i= 1, n1
!** accumulate squared differences
k(i,j)=  (x1(i,ic)- xtemp)**2 + k(i,j)
enddo
enddo
enddo
!**** at this point k( i,j) is the squared distance between x1_i and x2_j
!*** now evaluate radial basis functions
do  j =1,n2
do  i= 1, n1
k(i,j)=  radfun( k(i,j), par(1), par(2) )
enddo
enddo
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! evaluates thin plate spline radial basis function
double precision function radfun(d2, par1, par2)
double precision d2, par1, par2
if( d2.lt.1e-20) then
d2= 1e-20
endif
if( int(par2).eq.0) then
radfun= (d2)**( par1)
else
! note: d2 is squared distance
! divide by 2 to have log evaluated on distance
! as opposed to squared distance. 
radfun= (log(d2)/2) * ((d2)**( par1))
endif
return
end

! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html







subroutine rcss(h,npoint,x,y,wt,sy,trace,diag,cv, &
& ngrid,xg,yg,job,ideriv,din,dout,ierr)
! This a program to compute a robust univariate spline according to the
! model:
!   minimize   (1/n) sum(i=1,n)[rho( y(i)-f(x(i) )] + lambda*J(f)
!    over f
! definition of the rho function and its derivative are in rcsswt
! and rcssr
!
! One way of speeding convergence is to use the results from a previous
! estimate as the starting values for the another estimate. This is
! particulary appropriate when lambda or a parameter in the rho function
! is being varied. Moderate changes in lambda will often yield similar
! estimates. The way to take advantage of this is to pass the weights
! from the previous fit as teh starting values for the next estimate
!   
!   Arguments of rcss:   
!    h : natural log of lambda  
!                               
!        if h is passed with a value less than or equal -1000 no smoothing will   
!        be done and the spline will interploate the data points  
!    npoint: number of observations  
!    (x,y) : pairs of data points to be smoothed  
!            x(i) are assumed to be increasing. Repeated   
!            observations at the same x are not handled by this routine.  
!    sy : on return predicted values of f at x  
!    wt : weights used in the iterivatively reweighted least 
!           squares algorithm to compute robust spline. Vector
!           passed are used as the starting values. Subsequent
!           iterations compute the weights by a call to  the 
!           subroutine     rcsswt  
!
!          that is the linear approximation of teh estimator at 
!              convergence. 
!          trace= tr(A(lambda)) = " effective number of paramters"  
!   diag: diagonal elements of A(lambda) ( this is the most   
!         computationally intetnsive opertation in this subroutine.)  
!   cv: approximate cross-validation function 
!         using the linear approximation at convergence
! 
!   ngrid,xg,yg : on return, the ith deriv of the spline estimate is  
!                 evaluated on a grid, xg, with ngrid points.  The  
!                 values are returned in yg  
!  
!   ideriv: 0 = evaluate the function according to the job code  
!           1 = evaluate the first derivative of the function according  
!               to the job code  
!           2 = evaluate the second derivative of the function according  
!               to the job code  
! 
!   din:  Vector of input parameters 
!           din(1)= cost for cv
!           din(2)= offset for cv
!           din(3)= max number of iterations
!           din(4)= tolerance criterion for convergence
!
!           din(5)= C scale parameter in robust function (transition from
!                   quadratic to linear.
!           din(6)= alpha   1/2 slope of rho function for large, positive
!                   residuals   slope for residuals <0 : is  1/2 (1-alpha) 
!                  see comments in rcsswt for defintion of the  rho and psi 
!                  functions
!
!   job: in decimal job=(a,b,c)  (a=igcv, b=igrid)   
!        a=0  evaluate spline at x values, return predicted values in sy  
!        a=1  same as a=0 plus returning values of trace, diag and cv  
!        a=2  do no smoothing, interpolate the data  
!        a=3  same as a=1 but use the passed values in din(1) din(2)
!             for computing cv function
!  
!        b=0  do not evaluate the spline at any grid points  
!        b=1  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
!             of the spline at the (unique) sorted,data points, xg, return yg  
!        b=2  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
!             of the spline at ngrid equally spaced points between x(1)  
!             and x(npoints)      
!        b=3  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
!             of the spline at ngrid points on grid passed in xg.
!             NOTE: extrapolation of the estimate beyond the range of   
!                     the x's will just be a linear function.  
! 
!
!  
!         c=1    Assume that the X's are in sorted order
!         c=2 Do not sort the X's use the  current set of keys
!              Should only be used on a second call to smooth
!              same design

! Arguments of subroutine:
!    dout(1)= numit
!    dout(2)=tstop
!    dout(3) = trace
!    dout(4)= cv
!      numit: actual number of iterations for convergence
!      tstop: value of convergence criterion at finish.
!
! ierr: if ierr>0 something bad has happened
!       ierr>10 indicates problems in the call to the cubic spline
!       routine.
!  

parameter(NMAX=20000)
implicit double precision (a-h,o-z)
double precision h,trace,cv
double precision wt(npoint),X(npoint),Y(npoint)
double precision sy(npoint),diag(npoint)
double precision xg(ngrid),yg(ngrid)
double precision din(10), dout(10),cost,  offset, dum1, dum2

integer npoint,ngrid ,itj(3), job(3)


if(npoint.gt.NMAX) then
ierr=1
return
endif

maxit= int(din(3))
tstop=din(4)
ybar=0.0
ysd=0.0
do 10 k=1,npoint
diag(k) = y(k)
ybar= ybar+ y(k)
ysd= ysd+ y(k)**2
10 continue
ybar= ybar/npoint
ysd= sqrt( ysd/npoint - ybar**2)
! Start iterating
test=0.0

itj(1)= 0
itj(2)=0
itj(3)=0
do 500 it=1,maxit
if( it.gt.1) then
itj(3)=2
endif
! fit a weighted least squares cubic smoothing spline
call css(h,npoint,x,y,wt,sy, &
& dum1,diag,dum2,ngrid,xg,yg, &
& itj,ideriv,ierr)
! check for an error returned by spline routine
if(ierr.gt.0) then
!         add 10 so these can be distinguished from errors in this routine
ierr= 10 + ierr
return
endif

! compute convergence criterion
! The intent of this criterion is to find out when successive iterations
! produce changes that are small 


do 25 i=1,npoint
test=(diag(i)-sy(i))**2 + test
diag(i)= sy(i)
25 continue

test=sqrt(test/npoint)/ysd
if( test.lt.tstop  ) then
!            * exit loop *
numit=it
goto 1000
endif

!
!    make up new set of weights
!
call rcsswt( npoint, y,sy,wt, din(5))
!    reinitialize test criterion for convergence
!
test=0.0
500 continue

numit= maxit
1000 continue
! One last call if job code is not 0 
if( (job(1).ne.0).or.(job(2).ne.0)) then
!
call css(h,npoint,x,y,wt,sy, &
& trace,diag,cv,ngrid,xg,yg, &
& job,ideriv,ierr)

ia= job(1)
ig= job(2)
!      if(ig.gt.0)  then
!      endif
! calculate cv value if asked for 
if( (ia.eq.1) .or.( ia.eq.3) ) then
if(ia.eq.3) then
cost= din(1)
offset= din(2)/npoint
else
cost=1
offset= 0
endif

cv=0.0
do 1500 k=1,npoint
! compute approx. cross-validated residual


! plug cv residual into rho function, din(5) is the begining of parameter 
! vector for the rho function (scale and alpha)
!
! but only do this if the leverage is away from one. 
! this prevents the numerical problems with quantile splein of a zero
!  residual and 
! a zero denominator. 
if(1- diag(k).gt.1e-7) then
resid= (y(k)- sy(k))/( 1- cost*(diag(k)+offset))
cv= cv + rcssr(resid, din(5))
endif

1500 continue
cv= cv/npoint
endif
endif

dout(1)=numit
dout(2)=test
dout(3)=trace
dout(4)=cv

return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html


double precision function rcssr(r,par)
!
!     robust rho function:
!  This is a peicewise polynomial with knots at -C , 0 and C
!  the function is quadratic for -C<u<0 and 0<u<C
!  the function is linear for u<-C and u>C
!  rho is continuous for all u aqnd differentiable for all points 
!  except u=0 when a != 1/2 
!   
!
!    rho(u) =      2*a*u - a*c     for u>C
!                  a*u**2/C        for   0<u< C   
!                  (1-a)*u**2/C    for -C<u<0
!                  2*(1-a)*u - (1-a)*C  for u< -C
!
!        Note a= par(1), C= par(2)
implicit double precision (a-h, o-z)
double precision r, par(2),c,a
c= par(1)
if( r.gt.0 ) then
a=par(2)
else
a =(1-par(2))
r= -r
endif
if( r.le.c) then
rcssr= a*r*r/c
else
rcssr= 2*a*(r) - a*c
endif
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

!**********
subroutine rcsswt(n,y, sy, wt, par)
implicit double precision (a-h, o-z)
double precision y(n), sy(n), wt(n),psi,a,am1,c
double precision par(2)
!
!   psi(u) is the derivative of rho(u) defined in rcssr above
!
!   It is composed of peicewise linear and peicewise constant segements
!   and will be continuous except at u for a!=.5. 
!
!
a= par(2)
am1 = (1-par(2))
c= par(1)
do 100 k=1, n
!   find scaled residual
r= (y(k)- sy(k))/c
if( (r.gt. 0)) then
if( r.lt. 1) then
psi= 2*a*r

else
psi= 2*a
endif
else
if( r.gt.-1) then
psi= 2*am1*r

else
psi= -2*am1
endif

endif
!
! note weights supplied to cubic spline routine follow the convention that
! they are in terms of standard deviations. ( The more common form is
!   as reciprocal variances
!
wt(k) = dsqrt( 2*r/psi)
100 continue
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

subroutine rdist( nd,x1,n1,x2,n2, k)
integer nd,n1,n2,ic
double precision  x1(n1,nd), x2(n2,nd), k(n1,n2)
double precision xtemp

do  j =1,n2
xtemp= x2(j,1)
do   i= 1, n1
!** accumulate squared differences
k(i,j)=  (x1(i,1)- xtemp)**2
enddo
enddo
! **** loop through columns of output matrix K
!*** outer most loop over columns of x1 and x2 should reduce memory swaps
if( nd.ge.2) then
do  ic= 2, nd
do  j =1,n2
xtemp= x2(j,ic)
do   i= 1, n1
!** accumulate squared differences
k(i,j)=  (x1(i,ic)- xtemp)**2 + k(i,j)
enddo
enddo
enddo
endif
!**** at this point k( i,j) is the squared distance between x1_i and x2_j
do  j =1,n2
do   i= 1, n1
k(i,j)= sqrt( k(i,j))
enddo
enddo
return
end



subroutine rdist1( nd,x1,n1,k)
integer nd,n1,ic
double precision  x1(n1,nd), k(n1,n1)
double precision xtemp,  dtemp

do  j =1,n1
xtemp= x1(j,1)
do   i= 1, j
!** accumulate squared differences
k(i,j)=  (x1(i,1)- xtemp)**2
enddo
enddo
! **** loop through columns of output matrix K
!*** outer most loop over columns of x1 and x2 should reduce memory swaps
if( nd.ge.2) then
do  ic= 2, nd
do  j =1,n1
xtemp= x1(j,ic)
do   i= 1, j
!** accumulate squared differences
k(i,j)=  (x1(i,ic)- xtemp)**2 + k(i,j)
enddo
enddo
enddo
endif
!**** at this point k( i,j) is the squared distance between x1_i and x2_j
!**** for the upper triangle of matrix
do  j = 1,n1
do   i= 1, j
dtemp = sqrt( k(i,j))
k(i,j)= dtemp
!                
! filling in lower part takes a substantial time and is omitted  
! This means the returned matrix k has indeterminant vlaues in the
! lower triangle.              
!                k(j,i)= dtemp
enddo
enddo
return
end
! fields, Tools for spatial data
! Copyright (C) 2017, Institute for Mathematics Applied Geosciences
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html



! OLD routine with line numbers and computed if   
!       SUBROUTINE SORTM0(K,ki,N)  
!  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS K OF LENGTH N      
!  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.        
! integer array ki is permuted along with K  

!      double precision K(N),KK   
!      integer ki(N),kki  
!      INTEGER R               
!      IF(N.LE.1) RETURN      
!      L=N/2+1               
!      R=N                  
!  2   IF(L.GT.1) GO TO 1  
!      KK=K(R)  
!      kki= ki(R)                 
!      K(R)=K(1)  
!      ki(R)=ki(1)    
!      R=R-1         
!      IF(R.EQ.1) GO TO 9     
!      GO TO 3               
!  1   L=L-1                
!      KK=K(L)  
!      kki=ki(L)   
!  3   J=L         
!  4   I=J         
!      J=2*J       
!      IF(J-R) 5,6,8     
!  5   IF(K(J).LT.K(J+1)) J=J+1     
!  6   IF(KK.GT.K(J)) GO TO 8      
!  7   K(I)=K(J)  
!      ki(I)=ki(J)   
!      GO TO 4       
!  8   K(I)=KK      
!      ki(I)=kki   
!      GO TO 2    
!  9   K(1)=KK   
!      ki(1)=kki     
!      RETURN       
!      END

SUBROUTINE SORTM(X,ki,N)
!  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS X OF LENGTH N      
!  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.        
! integer array ki is permuted along with X    
double precision X(N),XX
integer N
integer ki(N),kki
INTEGER R, L,I, J
IF(N.LE.1) RETURN
L=N/2+1
R=N
2 IF(L.GT.1) GO TO 1
XX=X(R)
kki= ki(R)
X(R)=X(1)
ki(R)=ki(1)
R=R-1
IF(R.EQ.1) GO TO 9
GO TO 3
1 L=L-1
XX=X(L)
kki=ki(L)
3 J=L
4 I=J
J=2*J
!     IF(J-R) 5,6,8
! <  goes here      
if( (J-R).LT.0) then
IF(X(J).LT.X(J+1)) J=J+1
endif
! <= go to here       
if(  (J-R).LE.0) then
IF(XX.GT.X(J)) GO TO 8
X(I)=X(J)
ki(I)=ki(J)
GO TO 4
endif
! > goes to here       
8 X(I)=XX
ki(I)=kki
GO TO 2
9 X(1)=XX
ki(1)=kki
RETURN
END


end module fields_native
