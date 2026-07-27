! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

!------------------------------------------------- 
!     AUSP	---  (72) auspec, mulspe
!     COEFAB  ---  (72) fpec7 (74) canoca
!     CORNOM  ---  (72) autcor, mulcor, fftcor
!     CROSCO  ---  (72) autcor, mulcor
!     DMEADL  ---  (72) autcor, mulcor, fftcor
!     ECORCO  ---  (72) mulspe
!     FGER1	---  (72) mulnos, mulrsp
!     FGERCO  ---  (72) auspec, mulspe
!     INVDET  ---  (72) fpec7, optdes (74) canoca (78) mulmar, perars, mlomar, exsar
!     INVDETC ---  (72) mulnos, mulrsp
!     LTINV   ---  (72) wnoise (74) simcon
!     LTRVEC  ---  (72) wnoise (74) simcon
!     MATADL  ---  (72) fpec7, optdes
!     MIXRAD  ---  (72) fftcor (74) covgen
!     MULPLY  ---  (72) fpec7, optdes
!     MULVER  ---  (72) optsim (74) simcon
!     NEWSE   ---  (72) fpec7, (74) canoca
!     REARRA  ---  (72) fpec7
!     REARRAC ---  (72) mulfre
!     SIGNIF  ---  (72) auspec, mulspe
!     SMOSPE  ---  (72) mulspe
!     SUBD12  ---  (72) sglfre, mulfre
!     SUBDET  ---  (72) fpec7, (78) exsar
!     SUBNOS  ---  (72) mulnos
!     SUBTAC  ---  (72) fpec7, coefab (74) canoca
!     TRAMDL  ---  (72) fpec7, optdes
!     TRAMDR  ---  (72) fpec7, (74) canoca
!
!     DSUMF   ---  (72) autcor, mulcor, fftcor
!     RANDM   ---  (72) wnoise
!
!------------------------------------------------- 
!     INNERP  ---  (74) autarm, canoca, markov
!     MATINV  ---  (74) autarm, markov
!     MSVD    ---  (74) canarm, canoca
!     SUBTAC  ---  (74) canoca, simcon
!
!------------------------------------------------- 
!     ADDVAR  ---  (78) mulmar, perars, mlomar
!     AICCOM  ---  (78) mulmar, perars, mlomar
!     ARBAYS  ---  (78) unibar, blocar
!     ARCOEF  ---  (78) unibar, blocar, exsar
!     ARMFIT  ---  (78) unimar, bsubst, mlocar, exsar
!     BAYSPC  ---  (78) unibar, blocar
!     BAYSWT  ---  (78) unibar, mulbar, blocar, blomar
!     COEF2   ---  (78) mulmar, perars, mlomar
!     COMAIC  ---  (78) unimar, unibar, bsubst, mlocar, blocar, exsar
!     COPY    ---  (78) mulmar, mulbar, perars, mlocar, mlomar, blomar
!     DELETE  ---  (78) mulmar, perars, mlomar
!     FOUGER  ---  (78) unibar, mlocar, blocar, (72) raspec (74) nonst
!     HUSHL1  ---  (78) mulmar, mulmar, perars, mlomar, blomar
!     HUSHLD  ---  (78) unimar, unibar, bsubst, mulabr, perars, mlocar, blocar, mlomar, blomar, easer
!     MAICE   ---  (78) mulmar, perars, mlocar, blocar, mlomar, blomar, exsar
!     MARCOF  ---  (78) mulbar, blomar
!     MARFIT  ---  (78) mulmar, perars, mlomar
!     MBYSAR  ---  (78) mulbar, blomar
!     MBYSPC  ---  (78) mulbar, blomar
!     MCOEF   ---  (78) mulmar, perars, mlomar
!     MPARCO  ---  (78) mulbar, blomar
!     MRDATA  ---  (78) mulmar, mulbar, mlomar, blomar
!     MREDCT  ---  (78) mulmar, mulbar, perars, mlomar, blomar
!     MSDCOM  ---  (78) mulbar, blomar
!     MSETX1  ---  (78) mulmar, mulbar, perars, mlomar, blomar
!     NRASPE  ---  (78) unibar, bsubst, mlocar, blocar
!     PARCOR  ---  (78) exsar (84) decomp
!     RECOEF  ---  (78) unimar, bsubst, mlocar, exsar, (84) decomp
!     REDATA  ---  (78) unimar, unibar, bsubst, perars, mlocar, blocar, exsar
!     REDUCT  ---  (78) unimar, unibar, bsubst, mlocar, blocar, exsar
!     SDCOMP  ---  (78) unibar, bsubst, blocar
!     SETX1   ---  (78) unimar, unibar, bsubst, mlocar, blocar, exsar
!     SOLVE   ---  (78) mulbar, blomar
!     SRCOEF  ---  (78) mulmar, perars, mlomar
!     TRIINV  ---  (78) mulmar, perars
!
!     DMIN*8  ---  (78) blocar, blomar
!
!-------------------------------------------------
!
subroutine ausp(fc,p1,lagh1,a,la1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES SMOOTHED AUTO SPECTRUM.
!     FC: OUTPUT OF FGERCO
!     P1: SMOOTHED SPECTRUM
!     LAGH1: DIMENSION OF FC AND P1
!     A: SMOOTHING COEFFICIENTS
!     LA1: DIMENSION OF A (LESS THAN 11)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION FC(LAGH1),P1(LAGH1),A(LA1)
!      DIMENSION FC1(521)
!xx      DIMENSION FC1(LAGH1+2*(LA1-1))
integer lagh1, la1
real(dp) fc(lagh1), p1(lagh1), a(la1)
! local
integer la, lagshf
real(dp) fc1(lagh1+2*(la1-1))
!
la=la1-1
lagshf=lagh1+2*la
!     FC SHIFT-RIGHT BY LA FOR END CORRECTION
call ecorco(fc,lagh1,fc1,lagshf,la1)
!     SMOOTHING
call smospe(fc1,lagshf,a,la1,p1,lagh1)
return
end
!
!
!c	SUBROUTINE COEFAB(A1,B1,D,E,MS,K,MJ0,MJ)
subroutine coefab(a1,b1,d,e,ms,l,k)
  use timsac_kinds, only: dp
  implicit none
!     AR-FITTING
!     THIS SUBROUTINE COMPUTES FORWARD(A) AND BACKWARD(B) PREDICTOR
!     COEFFICIENTS.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	DIMENSION A1(MJ0,MJ,MJ),B1(MJ0,MJ,MJ)
!c	DIMENSION D(MJ,MJ),E(MJ,MJ)
!c	DIMENSION A(7,7),B(7,7),Z1(7,7),Z2(7,7)
!xx      DIMENSION A1(L,K,K),B1(L,K,K)
!xx      DIMENSION D(K,K),E(K,K)
!xx      DIMENSION A(K,K),B(K,K),Z1(K,K),Z2(K,K)
integer ms, l, k
real(dp) a1(l,k,k), b1(l,k,k), d(k,k), e(k,k)
! local
integer i, ii, jj, msm1, mmi
real(dp) a(k,k), b(k,k), z1(k,k), z2(k,k)
!
if(ms.eq.1) go to 40
msm1=ms-1
do 10 i=1,msm1
mmi=ms-i
!xx      DO 20 II=1,K
do 22 ii=1,k
do 20 jj=1,k
a(ii,jj)=a1(i,ii,jj)
!xx   20 B(II,JJ)=B1(MMI,II,JJ)
b(ii,jj)=b1(mmi,ii,jj)
20 continue
22 continue
!c	CALL MULPLY(D,B,Z1,K,K,K,MJ,MJ,MJ)
!c	CALL MULPLY(E,A,Z2,K,K,K,MJ,MJ,MJ)
!c	CALL SUBTAL(A,Z1,K,K,MJ,MJ)
!c	CALL SUBTAL(B,Z2,K,K,MJ,MJ)
call mulply(d,b,z1,k,k,k)
call mulply(e,a,z2,k,k,k)
call subtal(a,z1,k,k)
call subtal(b,z2,k,k)
!xx      DO 21 II=1,K
do 23 ii=1,k
do 21 jj=1,k
a1(i,ii,jj)=a(ii,jj)
!xx   21 B1(MMI,II,JJ)=B(II,JJ)
b1(mmi,ii,jj)=b(ii,jj)
21 continue
23 continue
10 continue
!xx   40 DO 30 II=1,K
40 do 31 ii=1,k
do 30 jj=1,k
a1(ms,ii,jj)=d(ii,jj)
!xx   30 B1(MS,II,JJ)=E(II,JJ)
b1(ms,ii,jj)=e(ii,jj)
30 continue
31 continue
return
end
!
subroutine cornom(c,cn,lagh1,cx0,cy0)
  use timsac_kinds, only: dp
  implicit none
!     NORMALIZATION OF COVARIANCE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION C(LAGH1),CN(LAGH1)
integer lagh1
real(dp) c(lagh1), cn(lagh1), cx0, cy0
! local
integer i
real(dp) cst1, ds
!
cst1=1.0d-00
ds=cst1/dsqrt(cx0*cy0)
do 10 i=1,lagh1
!xx   10 CN(I)=C(I)*DS
cn(i)=c(i)*ds
10 continue
return
end
!
subroutine crosco(x,y,n,c,lagh1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES C(L)=COVARIANCE(X(S+L),Y(S))
!     (L=0,1,...,LAGH1-1).
!xx      IMPLICIT REAL*8 (A-H,O-Z)
!xx      DIMENSION X(N),Y(N),C(LAGH1)
integer n, lagh1
real(dp) x(n), y(n), c(lagh1)
! local
integer i, ii, il, j, j1
real(dp) an, bn1, bn, ct0, t
!
an=n
bn1=1.0d-00
bn=bn1/an
ct0=0.0d-00
do 10 ii=1,lagh1
i=ii-1
t=ct0
il=n-i
do 20 j=1,il
j1=j+i
!xx   20 T=T+X(J1)*Y(J)
t=t+x(j1)*y(j)
20 continue
!xx   10 C(II)=T*BN
c(ii)=t*bn
10 continue
return
end
!
subroutine dmeadl(x,n,xmean)
  use timsac_kinds, only: dp
  implicit none
!     DOUBLE PRECISION MEAN DELETION
!c      DOUBLE PRECISION X,XMEAN,AN
!xx      DOUBLE PRECISION X,XMEAN,AN,DSUMF
!xx      DIMENSION X(N)
integer n
real(dp) x(n), xmean
! local
integer i
real(dp) an, dsumf
!
an=n
xmean=dsumf(x,n)/an
do 10 i=1,n
!xx   10 X(I)=X(I)-XMEAN
x(i)=x(i)-xmean
10 continue
return
end
!
subroutine ecorco(fc,lagh1,fc1,lagshf,la1)
  use timsac_kinds, only: dp
  implicit none
!     FC SHIFT-RIGHT BY LA FOR REAL PART END CORRECTION
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION FC(LAGH1),FC1(LAGSHF)
integer lagh1, lagshf, la1
real(dp) fc(lagh1), fc1(lagshf)
! local
integer i, i1, i2, i3, i4, la, la2, lagh2
!
lagh2=lagh1+1
la=la1-1
do 100 i=1,lagh1
i1=lagh2-i
i2=i1+la
!xx  100 FC1(I2)=FC(I1)
fc1(i2)=fc(i1)
100 continue
la2=lagh1+la
do 110 i=1,la
i1=la1-i
i2=la1+i
i3=la2-i
i4=la2+i
fc1(i1)=fc1(i2)
!xx  110 FC1(I4)=FC1(I3)
fc1(i4)=fc1(i3)
110 continue
return
end
!
!      SUBROUTINE FGER1
subroutine fger1(g,gr,gi,lg,h,jjf)
  use timsac_kinds, only: dp
  implicit none
!     FOURIER TRANSFORM(GOERTZEL METHOD)
!     THIS SUBROUTINE COMPUTES ONE VALUE OF THE FOURIER TRANSFORM BY
!     GOERTZEL METHOD.
!xx      IMPLICIT REAL*8(A-H,O-W)
!xx      INTEGER H
!      COMMON G,GR,GI,LG,H,JJF
!      DIMENSION G(31)
!xx      DIMENSION G(LG+1)
integer lg, h, jjf
real(dp) g(lg+1), gr, gi
! local
integer i, i2, lg3, lg4, lgp1
real(dp) cst0, t, pi, ah, ak, tk, ck, sk, ck2,&
&um2, um1, um0
cst0=0.0d-00
lgp1=lg+1
!     REVERSAL OF G(I),I=1,...,LGP1 INTO G(LG3-I)   LG3=LGP1+1
if(lgp1.le.1) go to 110
lg3=lgp1+1
lg4=lg3/2
do 100 i=1,lg4
i2=lg3-i
t=g(i)
g(i)=g(i2)
!xx  100 G(I2)=T
g(i2)=t
100 continue
110 pi=3.1415926536
ah=h
t=pi/ah
ak=jjf-1
tk=t*ak
ck=dcos(tk)
sk=dsin(tk)
ck2=ck+ck
um2=cst0
um1=cst0
if(lg.eq.0) go to 12
do 11 i=1,lg
um0=ck2*um1-um2+g(i)
um2=um1
!xx   11 UM1=UM0
um1=um0
11 continue
12 gr=ck*um1-um2+g(lgp1)
gi=-sk*um1
return
end
!
!
subroutine fgerco(g,lgp1,fc,lf1)
  use timsac_kinds, only: dp
  implicit none
!     FOURIER TRANSFORM (GOERTZEL METHOD)
!     THIS SUBROUTINE COMPUTES FOURIER TRANSFORM OF G(I),I=0,1,...,LG AT
!     FREQUENCIES K/(2*LF),K=0,1,...,LF AND RETURNS COSIN TRANSFORM IN
!     FC(K).
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION G(LGP1),FC(LF1)
integer lgp1, lf1
real(dp) g(lgp1), fc(lf1)
! local
integer i, i2, k, lg, lg3, lg4, lf
real(dp) t, pi, alf, ak, tk, ck, ck2,&
&um2, um1, um0
lg=lgp1-1
lf=lf1-1
!     REVERSAL OF G(I),I=1,...,LGP1 INTO G(LG3-I)   LG3=LGP1+1
if(lgp1.le.1) go to 110
lg3=lgp1+1
lg4=lgp1/2
do 100 i=1,lg4
i2=lg3-i
t=g(i)
g(i)=g(i2)
!xx  100 G(I2)=T
g(i2)=t
100 continue
110 pi=3.1415926536
alf=lf
t=pi/alf
do 10 k=1,lf1
ak=k-1
tk=t*ak
ck=dcos(tk)
ck2=ck+ck
um2=0.0d-00
um1=0.0d-00
if(lg.eq.0) go to 12
do 11 i=1,lg
um0=ck2*um1-um2+g(i)
um2=um1
!xx   11 UM1=UM0
um1=um0
11 continue
12 fc(k)=ck*um1-um2+g(lgp1)
10 continue
return
end
!
subroutine invdet(x,xdet,mm,mj)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!       THE INVERSE AND DETERMINANT OF X COMPUTATION                    
!                                                                       
!       INPUTS:                                                         
!          X:     MM*MM SQUARE MATRIX                                   
!          MM:    DIMENSION OF X                                        
!          MJ:    ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM           
!                                                                       
!       OUTPUTS:                                                        
!          X:     INVERSE OF X                                          
!          XDET:  DETERMINANT OF X                                      
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!xx      DIMENSION X(MJ,MJ)                                                
!c      DIMENSION  IDS(100)                                               
!xx      DIMENSION  IDS(MM)
integer mm, mj
real(dp) x(mj,mj), xdet
! local
integer i, j, jj, l, maxi, mm1, mmj, ids(mm)
real(dp) xmaxp, xc
!
xdet = 1.0d00
do 10 l=1,mm
!     PIVOTING AT L-TH STAGE                                            
xmaxp=0.10000d-10
maxi=0
do 110 i=l,mm
!xx    1 IF( DABS(XMAXP) .GE. DABS(X(I,L)) )     GO TO 110                 
if( dabs(xmaxp) .ge. dabs(x(i,l)) )     go to 110
xmaxp=x(i,l)
maxi=i
110 continue
ids(l)=maxi
if(maxi.eq.l) go to 120
if(maxi.gt.0) go to 121
xdet = 0.0d00
go to 140
!     ROW INTERCHANGE                                                   
121 do 14 j=1,mm
xc=x(maxi,j)
x(maxi,j)=x(l,j)
!xx   14 X(L,J)=XC
x(l,j)=xc
14 continue
xdet=-xdet
120 xdet=xdet*xmaxp
xc = 1.0d00 / xmaxp
x(l,l)=1.0d00
do 11 j=1,mm
!xx   11 X(L,J)=X(L,J)*XC
x(l,j)=x(l,j)*xc
11 continue
do 12 i=1,mm
if(i.eq.l) go to 12
xc=x(i,l)
x(i,l) = 0.0d00
do 13 j=1,mm
!xx   13 X(I,J)=X(I,J)-XC*X(L,J)                                           
x(i,j)=x(i,j)-xc*x(l,j)
13 continue
12 continue
10 continue
if(mm.gt.1) go to 123
go to 140
!     COLUMN INTERCHANGE                                                
123 mm1=mm-1
do 130 j=1,mm1
mmj=mm-j
jj=ids(mmj)
if(jj.eq.mmj) go to 130
do 131 i=1,mm
xc=x(i,jj)
x(i,jj)=x(i,mmj)
!xx  131 X(I,MMJ)=XC
x(i,mmj)=xc
131 continue
130 continue
140 return
end
!
!      SUBROUTINE INVDET(X,XDET,MM,MJ)
subroutine invdetc(x,xdet,mm)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES THE INVERSE AND DETERMINANT OF
!     UPPER LEFT MM X MM OF COMPLEX MATRIX X.
!     X: ORIGINAL MATRIX
!     MM: DIMENSION OF UPPER LEFT OF X (SHOULD BE LESS THAN 11)
!     XDET: DETERMINANT OF UPPER LEFT MM X MM OF X
!     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!     THE INVERSE MATRIX IS OVERWRITTEN ON THE ORIGINAL.
!      DIMENSION X(MJ,MJ)
!      DIMENSION IDS(10)
!xx      IMPLICIT COMPLEX*16(X)
!      DIMENSION X(MJ,MJ)
!      DIMENSION IDS(10)
!xx      DIMENSION X(MM,MM)
!xx      DIMENSION IDS(MM)
!xx      DOUBLE PRECISION CST0, CST1
integer mm
complex(kind(0.0d0)) x(mm,mm), xdet
! local
integer i, j, jj, l, mm1, mmj, maxi, ids(mm)
real(dp) cst0, cst1
complex(kind(0.0d0)) xmaxp, xc
!
cst0=0.0d-00
cst1=1.0d-00
xdet=cst1
do 10 l=1,mm
!     PIVOTING AT L-TH STAGE
xmaxp=0.10000d-10
maxi=0
do 110 i=l,mm
!xx      IF(CDABS(XMAXP).GE.CDABS(X(I,L))) GO TO 110
if(abs(xmaxp).ge.abs(x(i,l))) go to 110
xmaxp=x(i,l)
maxi=i
110 continue
ids(l)=maxi
if(maxi.eq.l) go to 120
if(maxi.gt.0) go to 121
xdet=cst0
go to 140
!     ROW INTERCHANGE
121 do 14 j=1,mm
xc=x(maxi,j)
x(maxi,j)=x(l,j)
!xx   14 X(L,J)=XC
x(l,j)=xc
14 continue
xdet=-xdet
120 xdet=xdet*xmaxp
xc=cst1/xmaxp
x(l,l)=cst1
do 11 j=1,mm
!xx   11 X(L,J)=X(L,J)*XC
x(l,j)=x(l,j)*xc
11 continue
do 12 i=1,mm
if(i.eq.l) go to 12
xc=x(i,l)
x(i,l)=cst0
do 13 j=1,mm
!xx   13 X(I,J)=X(I,J)-XC*X(L,J)
x(i,j)=x(i,j)-xc*x(l,j)
13 continue
12 continue
10 continue
if(mm.gt.1) go to 123
go to 140
!     COLUMN INTERCHANGE
123 mm1=mm-1
do 130 j=1,mm1
mmj=mm-j
jj=ids(mmj)
if(jj.eq.mmj) go to 130
do 131 i=1,mm
xc=x(i,jj)
x(i,jj)=x(i,mmj)
!xx  131 X(I,MMJ)=XC
x(i,mmj)=xc
131 continue
130 continue
140 return
end
!
!c      SUBROUTINE LTINV(R,K,MJ)                                          
subroutine ltinv(r,k)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE                                                 
!     THIS SUBROUTINE FACTORIZES (R(I,J): I,J=1,K) INTO R=L*L',         
!     WITH L LOWER TRIANGLE, AND GIVES L' ON AND ABOVE THE DIAGONAL OF R
!     MJ: ABSOLUTE DIMENSION OF R IN THE MAIN ROUTINE                   
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      DIMENSION R(MJ,MJ)                                                
!xx      DIMENSION R(K,K)
integer k
real(dp) r(k,k)
! local
integer i, l, l1, m
real(dp) cst1, rpivot, ril
!
cst1=1.0d-00
do 10 l=1,k
rpivot=cst1/dsqrt(r(l,l))
r(l,l)=cst1/rpivot
do 12 i=1,k
if(i.eq.l) go to 12
r(l,i)=rpivot*r(l,i)
12 continue
if(l.eq.k) go to 11
l1=l+1
do 13 i=l1,k
ril=-rpivot*r(i,l)
r(i,l)=ril*rpivot
do 14 m=1,k
if(m.eq.l) go to 14
r(i,m)=r(i,m)+ril*r(l,m)
14 continue
13 continue
10 continue
11 return
end
!
!c      SUBROUTINE LTRVEC(X,Y,Z,MM,NN,MJ1,MJ2)                            
subroutine ltrvec(x,y,z,mm,nn)
  use timsac_kinds, only: dp
  implicit none
!     Z=X*Y                                                             
!     (VECTOR Z)=(LOWER TRIANGLE OF UPPER LEFT MM X NN OF X)*(VECTOR Y) 
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE            
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      DIMENSION X(MJ1,MJ2),Y(NN),Z(MM)                                  
!xx      DIMENSION X(MM,NN),Y(NN),Z(MM)
integer mm, nn
real(dp) x(mm,nn), y(nn), z(mm)
! local
integer i, j
real(dp) cst0, sum
!
cst0=0.0d-00
do 10 i=1,mm
sum=cst0
do 11 j=1,i
!xx   11 SUM=SUM+X(I,J)*Y(J)
sum=sum+x(i,j)*y(j)
11 continue
!xx   10 Z(I)=SUM                                                          
z(i)=sum
10 continue
return
end
!
!      SUBROUTINE MATADL(X,Y,MM,NN,MJ1,MJ2)
subroutine matadl(x,y,mm,nn)
  use timsac_kinds, only: dp
  implicit none
!     MATRIX ADDITION
!     X=X+Y
!     (UPPER LEFT MM X NN OF X)=(UPPER LEFT MM X NN OF X)+(UPPER LEFT
!     MM X NN OF Y).
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X AND Y IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!      DIMENSION X(MJ1,MJ2),Y(MJ1,MJ2)
!xx      DIMENSION X(MM,NN),Y(MM,NN)
integer mm, nn
real(dp) x(mm,nn), y(mm,nn)
! local
integer i, j
!
!xx      DO 10 I=1,MM
do 11 i=1,mm
do 10 j=1,nn
!xx   10 X(I,J)=X(I,J)+Y(I,J)
x(i,j)=x(i,j)+y(i,j)
10 continue
11 continue
return
end
!
subroutine mixrad(z,n,n2p,isg)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE
!     MIXED RADIX FAST FOURIER TRANSFORM
!     ISG=-1...FOURIER TRANSFORM
!     ISG=1...INVERSE FOURIER TRANSFORM
!xx      IMPLICIT REAL*8(A-H,O-Y)
!xx      IMPLICIT COMPLEX*16(Z)
!xx      DIMENSION Z(N)
!c	DIMENSION MS(11)
!xx      DIMENSION MS(N2P)
integer n, n2p, isg
complex(kind(0.0d0)) z(n)
! local
integer i, i1, j, j1, j2, j3, j4, jm1, jf, l, ll, lm4, m, n3, n5,&
&m4, nm1, ms(n2p)
real(dp) cst0, cst1, an, pi, pi2, sg, am4, am5, ajm1, arg,&
&c1, s1, c2, s2, c3, s3
complex(kind(0.0d0)) zci, zw1, zw2, zw3, zc1, zc2, zc3, zc4, zc
cst0=0.0d-00
cst1=1.0d-00
an=n
pi=3.1415926536
pi2=pi+pi
sg=isg
!xx      ZCI=SG*DCMPLX(CST0,CST1)
zci=sg*cmplx(cst0,cst1,kind(0.0d0))
do 10 i=1,n2p
!xx   10 MS(I)=2**(N2P-I)
ms(i)=2**(n2p-i)
10 continue
n3=n2p/2
m=n
do 11 l=1,n3
m=m/4
m4=m*4
lm4=n-m4+1
am4=m4
am5=sg*pi2/am4
do 12 j=1,m
jm1=j-1
ajm1=jm1
arg=ajm1*am5
c1=dcos(arg)
s1=dsin(arg)
c2=c1*c1-s1*s1
s2=c1*s1+c1*s1
c3=c1*c2-s1*s2
s3=c1*s2+c2*s1
!xx      ZW1=DCMPLX(C1,S1)
!xx      ZW2=DCMPLX(C2,S2)
!xx      ZW3=DCMPLX(C3,S3)
zw1=cmplx(c1,s1,kind(0.0d0))
zw2=cmplx(c2,s2,kind(0.0d0))
zw3=cmplx(c3,s3,kind(0.0d0))
do 13 i=1,lm4,m4
j1=i+jm1
j2=j1+m
j3=j2+m
j4=j3+m
zc1=z(j1)+z(j3)
zc2=z(j1)-z(j3)
zc3=z(j2)+z(j4)
zc4=z(j2)-z(j4)
z(j1)=zc1+zc3
z(j2)=(zc1-zc3)*zw2
zc4=zci*zc4
z(j3)=(zc2+zc4)*zw1
z(j4)=(zc2-zc4)*zw3
13 continue
12 continue
11 continue
n5=n2p-2*n3
if(n5.ne.1) go to 120
nm1=n-1
do 110 i=1,nm1,2
i1=i+1
zc=z(i)+z(i1)
z(i1)=z(i)-z(i1)
z(i)=zc
110 continue
!     UNSCRAMBLING
120 jf=0
do 16 i=1,n
if(jf.lt.i) go to 17
zc=z(i)
z(i)=z(jf+1)
z(jf+1)=zc
17 do 18 l=1,n2p
ll=l
if(jf.lt.ms(l)) go to 19
!xx   18 JF=JF-MS(L)
jf=jf-ms(l)
18 continue
ll=n2p
19 jf=jf+ms(ll)
16 continue
if(isg.lt.0) go to 30
do 20 i=1,n
!xx   20 Z(I)=Z(I)/AN
z(i)=z(i)/an
20 continue
30 return
end
!
!      SUBROUTINE MULPLY(X,Y,Z,MM,NN,NC,MJ1,MJ2,MJ3)
subroutine mulply(x,y,z,mm,nn,nc)
  use timsac_kinds, only: dp
  implicit none
!     MATRIX MULTIPLICATION
!     Z=X*Y
!     (UPPER LEFT MM X NC OF Z)=(UPPER LEFT MM X NN OF X)*(UPPER LEFT
!     NN X NC OF Y).
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!     (MJ2,MJ3): ABSOLUTE DIMENSION OF Y IN THE MAIN ROUTINE
!     (MJ1,MJ3): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!      DIMENSION X(MJ1,MJ2),Y(MJ2,MJ3),Z(MJ1,MJ3)
!xx      DIMENSION X(MM,NN),Y(NN,NC),Z(MM,NC)
integer mm, nn, nc
real(dp) x(mm,nn), y(nn,nc), z(mm,nc)
! local
integer i, j, k
real(dp) cst0, sum
!
cst0=0.0d-00
do 10 i=1,mm
do 11 j=1,nc
sum=cst0
do 12 k=1,nn
!xx   12 SUM=SUM+X(I,K)*Y(K,J)
sum=sum+x(i,k)*y(k,j)
12 continue
z(i,j)=sum
11 continue
10 continue
return
end
!
!c      SUBROUTINE MULVER(X,Y,Z,MM,NN,MJ1,MJ2)                            
subroutine mulver(x,y,z,mm,nn)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE                                                 
!     Z=X*Y (X: MATRIX  Y,Z: VECTORS)                                   
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE            
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      DIMENSION X(MJ1,MJ2),Y(NN),Z(MM)                                  
!xx      DIMENSION X(MM,NN),Y(NN),Z(MM)
integer mm, nn
real(dp) x(mm,nn), y(nn), z(mm)
! local
integer i, j
real(dp) cst0, sum
!
cst0=0.0d-00
do 10 i=1,mm
sum=cst0
do 11 j=1,nn
!xx   11 SUM=SUM+X(I,J)*Y(J)                                               
sum=sum+x(i,j)*y(j)
11 continue
!xx   10 Z(I)=SUM                                                          
z(i)=sum
10 continue
return
end
!
!c	SUBROUTINE NEWSE(A1,SE,MS,K,MJ0,MJ)
subroutine newse(a1,cv,se,ms,l,k,lcv1)
  use timsac_kinds, only: dp
  implicit none
!     SE COMPUTATION
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	DIMENSION A1(MJ0,MJ,MJ)
!c	DIMENSION SE(MJ,MJ)
!c	DIMENSION CV(25,7,7)
!c	DIMENSION A(7,7),R(7,7),Z(7,7)
!xx      DIMENSION A1(L,K,K), CV(LCV1,K,K)
!xx      DIMENSION SE(K,K)
!xx      DIMENSION A(K,K),R(K,K),Z(K,K)
integer ms, l, k, lcv1
real(dp) a1(l,k,k), cv(lcv1,k,k), se(k,k)
! local
integer i, ii, jj, mmi, msp2
real(dp) a(k,k), r(k,k), z(k,k), cst0
!c	COMMON /COM10/CV
cst0=0.0d-00
!xx      DO 10 II=1,K
!xx      DO 10 JJ=1,K
!xx   10 Z(II,JJ)=CST0
z(1:k,1:k)=cst0
!
msp2=ms+2
do 11 i=1,ms
mmi=msp2-i
do 12 ii=1,k
!xx      DO 12 JJ=1,K
do 13 jj=1,k
a(ii,jj)=a1(i,ii,jj)
!xx   12 R(II,JJ)=CV(MMI,II,JJ)
r(ii,jj)=cv(mmi,ii,jj)
13 continue
12 continue
!c	CALL MULPLY(A,R,SE,K,K,K,MJ,MJ,MJ)
!c   11 CALL MATADL(Z,SE,K,K,MJ,MJ)
call mulply(a,r,se,k,k,k)
!xx   11 CALL MATADL(Z,SE,K,K)
call matadl(z,se,k,k)
11 continue
!
do 14 ii=1,k
!xx      DO 14 JJ=1,K
do 15 jj=1,k
!xx   14 R(II,JJ)=CV(MSP2,II,JJ)
r(ii,jj)=cv(msp2,ii,jj)
15 continue
14 continue
!c	CALL SUBTAC(R,Z,SE,K,K,MJ,MJ)
call subtac(r,z,se,k,k)
return
end
!
!      SUBROUTINE REARRA(X,INW,IP0,IP,MJ)
subroutine rearra(x,inw,ip0,ip)
  use timsac_kinds, only: dp
  implicit none
!     SUBMATRIX REARRANGEMENT
!     X: ORIGINAL MATRIX
!     INW: INDICATOR OF ADOPTED ROWS
!     IP0: DIMENSION OF ORIGINAL MATRIX, SHOULD BE LESS THAN 11
!     IP: DIMENSION OF REARRANGED SUBMATRIX
!     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!     THE REARRANGED SUBMATRIX IS OVERWRITTEN ON THE ORIGINAL.
!     NEXT STATEMENT SHOULD BE REPLACED BY
!     IMPLICIT COMPLEX*16(X)
!     FOR COMPLEX VERSION.
!xx      IMPLICIT REAL*8(X)
!      DIMENSION X(MJ,MJ),INW(IP)
!      DIMENSION IOD(10)
!xx      DIMENSION X(IP0,IP0),INW(IP)
!xx      DIMENSION IOD(IP0)
integer ip0, ip, inw(ip)
real(dp) x(ip0,ip0)
! local
integer i, i1, i2, ii, id, jj, iod(ip0)
real(dp) xc
!
do 300 i=1,ip0
!xx  300 IOD(I)=I
iod(i)=i
300 continue
do 301 i=1,ip
i1=inw(i)
i2=iod(i1)
if(i.eq.i2) go to 301
!     ROW INTERCHANGE
do 312 jj=1,ip0
xc=x(i,jj)
x(i,jj)=x(i2,jj)
!xx  312 X(I2,JJ)=XC
x(i2,jj)=xc
312 continue
!     COLUMN INTERCHANGE
do 314 ii=1,ip0
xc=x(ii,i)
x(ii,i)=x(ii,i2)
!xx  314 X(II,I2)=XC
x(ii,i2)=xc
314 continue
id=iod(i)
iod(i2)=id
iod(id)=i2
301 continue
return
end
!
subroutine rearrac(x,inw,ip0,ip)
  use timsac_kinds, only: dp
  implicit none
!     SUBMATRIX REARRANGEMENT
!     X: ORIGINAL MATRIX
!     INW: INDICATOR OF ADOPTED ROWS
!     IP0: DIMENSION OF ORIGINAL MATRIX, SHOULD BE LESS THAN 11
!     IP: DIMENSION OF REARRANGED SUBMATRIX
!     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!     THE REARRANGED SUBMATRIX IS OVERWRITTEN ON THE ORIGINAL.
!     NEXT STATEMENT SHOULD BE REPLACED BY
!xx      IMPLICIT COMPLEX*16(X)
!      DIMENSION X(MJ,MJ),INW(IP)
!      DIMENSION IOD(10)
!xx      DIMENSION X(IP0,IP0),INW(IP)
!xx      DIMENSION IOD(IP0)
integer ip0, ip, inw(ip)
complex(kind(0.0d0)) x(ip0,ip0)
! local
integer i, i1, i2, ii, id, jj, iod(ip0)
complex(kind(0.0d0)) xc
!
do 300 i=1,ip0
!xx  300 IOD(I)=I
iod(i)=i
300 continue
do 301 i=1,ip
i1=inw(i)
i2=iod(i1)
if(i.eq.i2) go to 301
!     ROW INTERCHANGE
do 312 jj=1,ip0
xc=x(i,jj)
x(i,jj)=x(i2,jj)
!xx  312 X(I2,JJ)=XC
x(i2,jj)=xc
312 continue
!     COLUMN INTERCHANGE
do 314 ii=1,ip0
xc=x(ii,i)
x(ii,i)=x(ii,i2)
!xx  314 X(II,I2)=XC
x(ii,i2)=xc
314 continue
id=iod(i)
iod(i2)=id
iod(id)=i2
301 continue
return
end
!
subroutine signif(p1,p2,p3,lagh1,n)
  use timsac_kinds, only: dp
  implicit none
!     SIGNIFICANCE TEST
!     P1: SPECTRUM SMOOTHED BY WINDOW W1
!     P2: SPECTRUM SMOOTHED BY WINDOW W2
!     P3: TEST STATISTICS
!     LAGH1: DIMENSION OF PI (I=1,2,3)
!     N: LENGTH OF THE ORIGINAL DATA
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION P1(LAGH1),P2(LAGH1),P3(LAGH1)
integer lagh1, n
real(dp) p1(lagh1), p2(lagh1), p3(lagh1)
! local
integer i, lagh
real(dp) h, an, han, sd2, sd3, t
lagh=lagh1-1
h=lagh
an=n
han=h/an
sd2=0.43d-00*dsqrt(han)
sd3=1.0d-00/sd2
do 10 i=1,lagh1
t=p2(i)/p1(i)-1.0d-00
!xx   10 P3(I)=DABS(T)*SD3
p3(i)=dabs(t)*sd3
10 continue
return
end
!
subroutine smospe(x,lagshf,a,la1,z,lagh1)
  use timsac_kinds, only: dp
  implicit none
!     SPECTRUM SMOOTHING BY THE FORMULA
!     Z(I)=A(0)X(I)+A(1)(X(I+1)+X(I-1))+...+A(LA)(X(I+LA)+X(I-LA))
!     I=0,1,...,LAGH.
!     ACTUAL X(I) IS SHIFTED TO THE RIGHT BY LA FOR END CORRECTION.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION X(LAGSHF),A(LA1),Z(LAGH1)
integer lagshf, la1, lagh1
real(dp) x(lagshf), a(la1), z(lagh1)
! local
integer i, i0, j, j1, j2, la
real(dp) sum1
!
la=la1-1
do 10 i=1,lagh1
i0=i+la
sum1=0.0d-00
do 11 j=1,la
j1=i0-j
j2=i0+j
!xx   11 SUM1=SUM1+A(J+1)*(X(J1)+X(J2))
sum1=sum1+a(j+1)*(x(j1)+x(j2))
11 continue
!xx   10 Z(I)=A(1)*X(I0)+SUM1
z(i)=a(1)*x(i0)+sum1
10 continue
return
end
!
subroutine subd12(n,lagh,k,d1,d2)
  use timsac_kinds, only: dp
  implicit none
!     CONSTANTS D1,D2 COMPUTATION
!xx      IMPLICIT REAL*8(A-H,O-Z)
!     L1: NUMBER OF A(I)S (LESS THAN 5)
!xx      DIMENSION A(4)
integer n, lagh, k
real(dp) d1, d2
! local
integer i, l1, nf
real(dp) a(4), cst0, an, h, sum, fk, c1, c2
cst0=0.0d-00
l1=2
a(1)=0.5d-00
a(2)=0.25d-00
an=n
h=lagh
sum=0.0d-00
do 20 i=2,l1
!xx   20 SUM=SUM+A(I)**2
sum=sum+a(i)**2
20 continue
sum=sum+sum+a(1)**2
sum=sum+sum
!xx      NF=AN/(H*SUM)+0.5D-00
nf=int(an/(h*sum)+0.5d-00)
fk=nf-k
if(fk.eq.cst0) go to 100
c1=fk-1.40d-00
if(c1.eq.cst0) go to 100
d1=(3.84d-00+10.0d-00/c1)/fk
if(d1.lt.cst0) go to 100
d1=dsqrt(d1)
go to 110
100 d1=100.0d-00
110 c2=fk+fk-1.40d-00
if(c2.eq.cst0) go to 120
d2=(3.0d-00+10.0d-00/c2)/fk
if(d2.lt.cst0) go to 120
d2=dsqrt(d2)
go to 130
120 d2=100.0d-00
130 return
end
!
subroutine subdet(x,xdetmi,mm,mj)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!       DETERMINANT OF X COMPUTATION                                    
!                                                                       
!     THIS SUBROUTINE COMPUTES THE DETERMINANT OF UPPER LEFT MM X MM    
!     OF X.  FOR GENERAL USE STATEMENTS 20-21 SHOULD BE RESTORED.       
!     X: ORIGINAL MATRIX                                                
!     XDETMI: DETERMINANT OF UPPER LEFT MM X MM OF X                    
!     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE                   
!xx      IMPLICIT REAL*8(X)                                                
!xx      DIMENSION X(MJ,MJ)
integer mm, mj
real(dp) x(mj,mj), xdetmi
! local
integer i, i1, j, jj, k, mm1
real(dp) xc, xxc
xdetmi=1.0d0
if(mm.eq.1) go to 18
mm1=mm-1
do 10 i=1,mm1
!xx   20 IF(X(I,I).NE.0.0D0) GO TO 11                                      
if(x(i,i).ne.0.0d0) go to 11
do 12 j=i,mm
if(x(i,j).eq.0.0d0) go to 12
jj=j
go to 13
12 continue
xdetmi=0.0d0
go to 17
13 do 14 k=i,mm
xxc=x(k,jj)
x(k,jj)=x(k,i)
!xx   14 X(K,I)=XXC
x(k,i)=xxc
14 continue
!xx   21 XDETMI=-XDETMI
xdetmi=-xdetmi
11 xdetmi=xdetmi*x(i,i)
xc=1.0d0/x(i,i)
i1=i+1
do 15 j=i1,mm
xxc=x(j,i)*xc
do 16 k=i1,mm
!xx   16 X(J,K)=X(J,K)-X(I,K)*XXC
x(j,k)=x(j,k)-x(i,k)*xxc
16 continue
15 continue
10 continue
18 xdetmi=xdetmi*x(mm,mm)
17 return
end
!
subroutine subnos(x,sd,ip,rs,r,mj)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES RELATIVE POWER CONTRIBUTIONS.
!     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!     IP: DIMENSION OF RS1 OR RL (LESS THAN 11)
!xx      IMPLICIT REAL*8(A-H,O-W)
!xx      IMPLICIT COMPLEX*16(X-Z)
!xx      DIMENSION X(MJ,MJ)
!xx      DIMENSION SD(MJ,MJ),RS(MJ,MJ),R(MJ,MJ)
!      DIMENSION RS1(10),RL(10)
!xx      DIMENSION RS1(MJ),RL(MJ)
integer ip, mj
real(dp) sd(mj,mj), rs(mj,mj), r(mj,mj)
complex(kind(0.0d0)) x(mj,mj)
! local
integer ii, jj, ll
real(dp) rs1(mj), rl(mj), cst0, cst1, sum, rx, rix, rconst
!
cst0=0.0d-00
cst1=1.0d-00
do 10 ii=1,ip
sum=cst0
do 11 jj=1,ip
!xx      RX=DREAL(X(II,JJ))
!xx      RIX=DIMAG(X(II,JJ))
rx=real(x(ii,jj))
rix=aimag(x(ii,jj))
rs1(jj)=(rx**2+rix**2)*sd(jj,jj)
sum=sum+rs1(jj)
!xx   11 RL(JJ)=SUM
rl(jj)=sum
11 continue
rconst=cst1/rl(ip)
do 14 jj=1,ip
!xx   14 RS(II,JJ)=RS1(JJ)*RCONST
rs(ii,jj)=rs1(jj)*rconst
14 continue
do 12 ll=1,ip
!xx   12 R(II,LL)=RL(LL)*RCONST
r(ii,ll)=rl(ll)*rconst
12 continue
10 continue
return
end
!
!c	SUBROUTINE SUBTAL(X,Y,MM,NN,MJ1,MJ2)
subroutine subtal(x,y,mm,nn)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE
!     MATRIX SUBTRACTION
!     X=X-Y
!     (UPPER LEFT MM X NN OF X)=(UPPER LEFT MM X NN OF X)-(UPPER LEFT
!     MM X NN OF Y).
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X AND Y IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	DIMENSION X(MJ1,MJ2),Y(MJ1,MJ2)
!xx      DIMENSION X(MM,NN),Y(MM,NN)
integer mm, nn
real(dp) x(mm,nn), y(mm,nn)
! local
integer i, j
!
!xx      DO 10 I=1,MM
do 11 i=1,mm
do 10 j=1,nn
!xx   10 X(I,J)=X(I,J)-Y(I,J)
x(i,j)=x(i,j)-y(i,j)
10 continue
11 continue
return
end
!
!      SUBROUTINE TRAMDL(X,Y,Z,MM,NN,NC,MJ1,MJ2,MJ3)
subroutine tramdl(x,y,z,mm,nn,nc)
  use timsac_kinds, only: dp
  implicit none
!     TRANSPOSE MULTIPLY (LEFT)
!     Z=X'*Y
!     (UPPER LEFT NN X NC OF Z)=(UPPER LEFT MM X NN OF X)'*(UPPER LEFT
!     MM X NC OF Y).
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!     (MJ1,MJ3): ABSOLUTE DIMENSION OF Y IN THE MAIN ROUTINE
!     (MJ2,MJ3): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!      DIMENSION X(MJ1,MJ2),Y(MJ1,MJ3),Z(MJ2,MJ3)
!xx      DIMENSION X(MM,NN),Y(MM,NC),Z(NN,NC)
integer mm, nn, nc
real(dp) x(mm,nn), y(mm,nc), z(nn,nc)
! local
integer i, j, k
real(dp) cst0, sum
!
cst0=0.0d-00
do 10 i=1,nn
do 11 j=1,nc
sum=cst0
do 12 k=1,mm
!xx   12 SUM=SUM+X(K,I)*Y(K,J)
sum=sum+x(k,i)*y(k,j)
12 continue
z(i,j)=sum
11 continue
10 continue
return
end
!
!c	SUBROUTINE TRAMDR(X,Y,Z,MM,NN,NC,MJ1,MJ2,MJ3)
subroutine tramdr(x,y,z,mm,nn,nc)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE
!     TRANSPOSE MULTIPLY (RIGHT)
!     Z=X*Y'
!     (UPPER LEFT MM X NC OF Z)=(UPPER LEFT MM X NN OF X)*(UPPER LEFT
!     NC X NN OF Y)'.
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!     (MJ3,MJ2): ABSOLUTE DIMENSION OF Y IN THE MAIN ROUTINE
!     (MJ1,MJ3): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	DIMENSION X(MJ1,MJ2),Y(MJ3,MJ2),Z(MJ1,MJ3)
!xx      DIMENSION X(MM,NN),Y(NC,NN),Z(MM,NC)
integer mm, nn, nc
real(dp) x(mm,nn), y(nc,nn), z(mm,nc)
! local
integer i, j, k
real(dp) cst0, sum
!
cst0=0.0d-00
do 10 i=1,mm
do 11 j=1,nc
sum=cst0
do 12 k=1,nn
!xx   12 SUM=SUM+X(I,K)*Y(J,K)
sum=sum+x(i,k)*y(j,k)
12 continue
z(i,j)=sum
11 continue
10 continue
return
end
!
!
real(dp) function dsumf(x,n)
  use timsac_kinds, only: dp
  implicit none
!     DOUBLE PRECISION SUMMATION
!xx      DOUBLE PRECISION X
!xx      DIMENSION X(N)
integer n
real(dp) x(n)
integer i
dsumf=0.0d-00
do 10 i=1,n
!xx   10 DSUMF=DSUMF+X(I)
dsumf=dsumf+x(i)
10 continue
return
end
!
!c	FUNCTION RANDOM(K)
!cc	 FUNCTION RANDM(K)
real(dp) function randm(k,k1,k2,k3,k4)
  use timsac_kinds, only: dp
  implicit none
!     RANDOM NUMBER GENERATOR
integer k, k1, k2, k3, k4
! local
integer j, m1, m2, m3, m4, mcst11, mc100
real(dp) x1, x2
!
mcst11=11
mc100=100
!c      IF(K) 1,2,1
if(k .lt. 0) go to 1
if(k .eq. 0) go to 2
if(k .gt. 0) go to 1
!     STARTING NUMBER FOR GENERATOR
1 k1=53
k2=95
k3=27
k4=04
!c	WRITE(6,4) K1,K2,K3,K4
2 m1=mcst11*k4
m2=mcst11*k3
m3=mcst11*k2+k4
m4=mcst11*k1+k3
j=m1/mc100
k4=m1-mc100*j
m2=m2+j
j=m2/mc100
k3=m2-mc100*j
m3=m3+j
j=m3/mc100
k2=m3-mc100*j
m4=m4+j
j=m4/mc100
k1=m4-mc100*j
x1=k1
x2=k2
!c	RANDOM=X1*1.E-2+X2*1.E-4
!cc      RANDM=X1*1.E-2+X2*1.E-4
randm=x1*1.d-2+x2*1.d-4
return
!xx    4 FORMAT(1H ,29HSTARTING NUMBER FOR GENERATOR,5X,4I3)
end
!
!
!------------------------------------------------- TIMSAC74
subroutine innerp(dd1,dd2,dinp12,inp)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE
!     INNER-PRODUCT OF DD1 AND DD2.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION DD1(INP),DD2(INP)
integer inp
real(dp) dd1(inp), dd2(inp), dinp12
! local
integer i
real(dp) cst0, sum
!
cst0=0.0d-00
sum=cst0
do 100 i=1,inp
!xx  100 SUM=SUM+DD1(I)*DD2(I)
sum=sum+dd1(i)*dd2(i)
100 continue
dinp12=sum
return
end
!
!
!c	SUBROUTINE MATINV(DET,M)
!xx      SUBROUTINE MATINV(M,HS,NN,IFG,LU)
subroutine matinv(m,hs,nn)
  use timsac_kinds, only: dp
  implicit none
!      HS IS AN M*M MATRIX (IN COMMON AREA).
!     HS-INVERSE IS RETURNED IN HS.
!     DETERMINANT IS RETURNED IN DET.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	COMMON /COM50/HS
!c	DIMENSION HS(50,50)
!xx      DIMENSION HS(NN,NN)
integer m, nn
real(dp) hs(nn,nn), xdet
!
call invdet(hs,xdet,m,nn)
!
!c      CST0=0.0D-00
!c      CST1=1.0D-00
!c      DET=CST1
!c      DO 1 J=1,M
!c      PVT=HS(J,J)
!c      DET=DET*PVT
!c	WRITE(6,60) J,PVT,DET
!c      IF (IFG.NE.0) WRITE(LU,60) J,PVT,DET
!c      IF (PVT.GT.0) GO TO 5
!c	WRITE(6,61)
!c      IF (IFG.NE.0) WRITE(LU,61)
!c    5 HS(J,J)=CST1
!c      DO 2 K=1,M
!c    2 HS(J,K)=HS(J,K)/PVT
!c      DO 1 K=1,M
!c      IF(K-J) 3,1,3
!c    3 T=HS(K,J)
!c      HS(K,J)=CST0
!c      DO 4 L=1,M
!c    4 HS(K,L)=HS(K,L)-HS(J,L)*T
!c    1 CONTINUE
!c   60 FORMAT(1H ,'J=',I5,5X,'PVT=',D12.5,5X,'DET=',D12.5)
!c   61 FORMAT(1H ,'WARNING: NON POSITIVE DEFINITE HESSIAN')
return
end
!
subroutine msvd(u,v,q,m,n,mj2,mj1)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE
!     THIS SUBROUTINE COMPLETES THE SINGULAR VALUE
!     DECOMPOSITION OF A REAL RECTANGULAR MATRIX A INTO THE FORM
!     A=U*DIAG(Q)*V' WITH ORTHGONAL MATRICES U AND V.
!     INPUTS:
!     U: ORIGINAL MATRIX A
!     M: NUMBER OF ROWS OF A, NOT LESS THAN N.
!     N: NUMBER OF COLUMNS OF A, NOT GREATER THAN 104.
!     OUTPUTS:
!     V: ORTHOGONAL MATRIX V
!     Q: SINGULAR VALUES IN DECREASING ORDER
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	DIMENSION E(50)
!xx      DIMENSION U(MJ2,MJ1),V(MJ1,MJ1)
!xx      DIMENSION E(M)
!xx      DIMENSION Q(N)
integer m, n, mj2, mj1
real(dp) u(mj2,mj1), v(mj1,mj1), q(n)
! local
integer i, ii, ip1, j, k, kk, kp1, l, ll, lp1, l1, nm1, np1
real(dp) e(m), eps, tol, cst0, cst1, cst2, c, f, g, fg, h,&
&s, t, x, y, z
!     ESP: SMALL POSITIVE CONSTANT TO BE USED FOR THE DECISION OF CONVER
!     TOL: (SMALLEST POSITIVE NUMBER REPRESENTABLE IN THE COMPUTER)/EPS
eps=1.0d-15
tol=1.0d-60
!     HOUSEHOLDER'S REDUCTION TO BIDIAGONAL FORM
cst0=0.0d-00
cst1=1.0d-00
cst2=2.0d-00
g=cst0
x=cst0
do 10 i=1,n
e(i)=g
s=cst0
l=i+1
do 11 j=i,m
!xx   11 S=S+U(J,I)**2
s=s+u(j,i)**2
11 continue
if(s.ge.tol) go to 12
g=cst0
go to 19
12 f=u(i,i)
!xx   13 IF(F.LT.CST0) GO TO 14
if(f.lt.cst0) go to 14
g=-dsqrt(s)
go to 15
14 g=dsqrt(s)
15 h=f*g-s
u(i,i)=f-g
if(l.gt.n) go to 19
do 16 j=l,n
s=cst0
do 17 k=i,m
!xx   17 S=S+U(K,I)*U(K,J)
s=s+u(k,i)*u(k,j)
17 continue
f=s/h
do 18 k=i,m
!xx   18 U(K,J)=U(K,J)+F*U(K,I)
u(k,j)=u(k,j)+f*u(k,i)
18 continue
16 continue
19 q(i)=g
s=cst0
if(l.gt.n) go to 201
do 20 j=l,n
!xx   20 S=S+U(I,J)**2
s=s+u(i,j)**2
20 continue
201 if(s.ge.tol) go to 22
g=cst0
go to 30
22 f=u(i,i+1)
!xx   23 IF(F.LT.CST0) GO TO 24
if(f.lt.cst0) go to 24
g=-dsqrt(s)
go to 25
24 g=dsqrt(s)
25 h=f*g-s
u(i,i+1)=f-g
if(l.gt.n) go to 30
do 26 j=l,n
!xx   26 E(J)=U(I,J)/H
e(j)=u(i,j)/h
26 continue
do 27 j=l,m
s=cst0
do 28 k=l,n
!xx   28 S=S+U(J,K)*U(I,K)
s=s+u(j,k)*u(i,k)
28 continue
do 29 k=l,n
!c   29 U(J,K)=U(J,K)+S*E(K)
u(j,k)=u(j,k)+s*e(k)
29 continue
27 continue
30 y=dabs(q(i))+dabs(e(i))
!     X=MAX.Y
if(y.gt.x) x=y
10 continue
np1=n+1
l=n+1
!     ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS
do 110 ii=1,n
i=np1-ii
if(l.gt.n) go to 202
if(g.eq.cst0) go to 115
h=u(i,i+1)*g
do 111 j=l,n
!xx  111 V(J,I)=U(I,J)/H
v(j,i)=u(i,j)/h
111 continue
do 112 j=l,n
s=cst0
do 113 k=l,n
!xx  113 S=S+U(I,K)*V(K,J)
s=s+u(i,k)*v(k,j)
113 continue
do 114 k=l,n
!xx  114 V(K,J)=V(K,J)+S*V(K,I)
v(k,j)=v(k,j)+s*v(k,i)
114 continue
112 continue
115 do 116 j=l,n
v(i,j)=cst0
!xx  116 V(J,I)=CST0
v(j,i)=cst0
116 continue
202 v(i,i)=cst1
g=e(i)
l=i
110 continue
!     DIAGONALIZATION OF THE BIDIAGONAL FORM
!xx  129 EPS=EPS*X
eps=eps*x
do 130 kk=1,n
k=np1-kk
kp1=k+1
!     TEST F SPLITTING
145 do 131 ll=1,k
l=kp1-ll
if(dabs(e(l)).le.eps) go to 134
if(dabs(q(l-1)).le.eps) go to 132
131 continue
!     CANCELLATION OF E(L) IF L>1.
132 c=cst0
s=cst1
l1=l-1
do 133 i=l,k
f=s*e(i)
e(i)=c*e(i)
if(dabs(f).le.eps) go to 134
g=q(i)
h=dsqrt(f*f+g*g)
q(i)=h
c=g/h
s=-f/h
133 continue
!     TEST F CONVERGENCE
134 z=q(k)
if(l.eq.k) go to 142
!     SHIFT FROM BOTTOM 2X2 MINOR
x=q(l)
y=q(k-1)
g=e(k-1)
h=e(k)
f=((y-z)*(y+z)+(g-h)*(g+h))/(cst2*h*y)
t=f*f+cst1
g=dsqrt(t)
if(f.ge.cst0) go to 135
fg=f-g
go to 136
135 fg=f+g
136 f=((x-z)*(x+z)+h*(y/fg-h))/x
!     QR TRANSFORMATION
!     NEXT QR TRANSFORMATION
c=cst1
s=cst1
lp1=l+1
do 137 i=lp1,k
g=e(i)
y=q(i)
h=s*g
g=c*g
z=dsqrt(f*f+h*h)
e(i-1)=z
c=f/z
s=h/z
f=x*c+g*s
g=-x*s+g*c
h=y*s
y=y*c
do 138 j=1,n
x=v(j,i-1)
z=v(j,i)
v(j,i-1)=x*c+z*s
!xx  138 V(J,I)=-X*S+Z*C
v(j,i)=-x*s+z*c
138 continue
!xx  139 Z=DSQRT(F*F+H*H)
z=dsqrt(f*f+h*h)
q(i-1)=z
c=f/z
s=h/z
f=c*g+s*y
x=-s*g+c*y
137 continue
e(l)=cst0
e(k)=f
q(k)=x
go to 145
!     CONVERGENCE
142 if(z.ge.cst0) go to 130
!     Q(K) IS MADE NON-NEGATIVE.
q(k)=-z
do 144 j=1,n
!xx  144 V(J,K)=-V(J,K)
v(j,k)=-v(j,k)
144 continue
130 continue
!     Q AND V ARE ORDERED IN DECREASING ORDER OF THE SINGULAR VALUES.
if(n.le.1) go to 340
nm1=n-1
do 311 i=1,nm1
ip1=i+1
do 312 j=ip1,n
if(q(i).ge.q(j)) go to 312
t=q(i)
q(i)=q(j)
q(j)=t
do 321 l=1,n
t=v(l,i)
v(l,i)=v(l,j)
v(l,j)=t
321 continue
312 continue
311 continue
340 return
end
!

!
!c	SUBROUTINE SUBTAC(X,Y,Z,MM,NN,MJ1,MJ2)
subroutine subtac(x,y,z,mm,nn)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE
!     MATRIX SUBTRACTION
!     Z=X-Y
!     (UPPER LEFT MM X NN OF Z)=(UPPER LEFT MM X NN OF X)-(UPPER LEFT
!     MM X NN OF Y).
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X, Y AND Z IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	DIMENSION X(MJ1,MJ2),Y(MJ1,MJ2),Z(MJ1,MJ2)
!xx      DIMENSION X(MM,NN),Y(MM,NN),Z(MM,NN)
integer mm, nn
real(dp) x(mm,nn), y(mm,nn), z(mm,nn)
! local
integer i, j
!xx      DO 10 I=1,MM
do 20 i=1,mm
do 10 j=1,nn
!xx   10 Z(I,J)=X(I,J)-Y(I,J)
z(i,j)=x(i,j)-y(i,j)
10 continue
20 continue
return
end
!
!-------------------------------------------------------------------
!
!c      SUBROUTINE  ADDVAR( X,D,IND,JND,K,L,M,MJ )                        
subroutine  addvar( x,ind,jnd,k,l,m,mj )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!         +-----------------------------------------------------------+ 
!         ! ADDITION OF THE VARIABLE M AS THE JJ-TH REGRESSOR (JJ<=L) ! 
!         +-----------------------------------------------------------+ 
!                                                                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
!             HUSHL1                                                    
!       ----------------------------------------------------------------
!       INPUTS:                                                         
!          X:            (K+1)*(K+1)  MATRIX                            
!          D:            WORKING AREA                                   
!          IND:          WORKING AREA                                   
!          JND(I)=J:     PRESENT STATUS, I-TH REGRESSOR IS VARIABLE J.  
!          K:            NUMBER OF VARIABLES                            
!          L:            NUMBER OF REGRESSORS IN THE PRESENT MODEL      
!          M:            INDICATION OF THE VARIABLE TO BE ADDED         
!          MJ:           ABSOLUTE DIMENSION OF X                        
!                                                                       
!       OUTPUTS:                                                        
!          X:            (K+1)*(K+1)  MATRIX                            
!          JND(I)=J:     UPDATED STATUS, I-TH REGRESSOR IS VARIABLE J.  
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  X(MJ,1) , D(1) , IND(1) , JND(1)                       
!x      DIMENSION  X(MJ,1) , IND(1) , JND(1)                       
!xx      DIMENSION  X(MJ,K+1) , IND(K+1) , JND(K+1)                       
integer k, l, m, mj, ind(k+1), jnd(k+1)
real(dp) x(mj,k+1)
! local
integer i, i1, ii1, j, jj, k1
!                                                                       
k1 = k + 1
do  60     i=1,k1
j = jnd(i)
!xx   60 IND(J) = I                                                        
ind(j) = i
60 continue
jj = ind(m)
!c      IF( JJ-L )     40,40,10                                           
if( jj-l .lt. 0 )  go to 40
if( jj-l .eq. 0 )  go to 40
if( jj-l .gt. 0 )  go to 10
10 ii1 = l + 1
do  20     i=ii1,jj
i1 = jj+l-i
!xx   20 JND(I1+1) = JND(I1)                                               
jnd(i1+1) = jnd(i1)
20 continue
jnd(l) = m
!                                                                       
!c      CALL  HUSHL1( X,D,MJ,K1,JJ,L,IND,JND )                            
call  hushl1( x,mj,k1,jj,l,ind,jnd )
!                                                                       
!xx   30 L = L + 1                                                         
l = l + 1
40 return
!                                                                       
end
!
!
subroutine  aiccom( x,n,m,k,mj,sd,aic )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE COMPUTES INNOVATION VARIANCE AND AIC OF THE MODEL 
!     WITH M REGRESSORS.                                                
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  X(MJ,1)                                                
!xx      DIMENSION  X(MJ,K+1)
integer n, m, k, mj
real(dp) x(mj,k+1), sd, aic
! local
integer i, k1, m1
real(dp) fn, sum
!                                                                       
!       INPUTS:                                                         
!          X:     (K+1)*(K+1) REDUCED MATRIX                            
!          N:     DATA LENGTH                                           
!          M:     NUMBER OF REGRESSORS OF THE MODEL                     
!          K:     POSSIBLE LARGEST NUMBER OF REGRESSORS                 
!          MJ:    ABSOLUTE DIMENSION OF X                               
!                                                                       
!       OUTPUTS:                                                        
!          SD:    INNOVATION VARIANCE OF THE MODEL                      
!          AIC:   AIC OF THE MODEL                                      
!                                                                       
fn = n
m1 = m + 1
k1 = k + 1
sum = 0.d0
do 10  i=m1,k1
!xx   10 SUM = SUM + X(I,K1)*X(I,K1)                                       
sum = sum + x(i,k1)*x(i,k1)
10 continue
sd = sum / fn
!xx   20 AIC = FN*DLOG( SD ) + 2.D0*M                                      
aic = fn*dlog( sd ) + 2.d0*m
return
end
!
!
!c      SUBROUTINE  ARBAYS( X,D,K,LAG,N,ISW,TITLE,MJ1,A,B,SDB,AICB )      
!xx      SUBROUTINE  ARBAYS( X,D,K,LAG,N,ISW,MJ1,SD,AIC,DIC,AICM,SDMIN,
subroutine  arbays( x,d,k,n,isw,mj1,sd,aic,dic,aicm,sdmin,&
&imin,a,b1,b,c,sdb,pn,aicb )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!         +-----------------------------------------------+             
!         ! AUTOREGRESSIVE MODEL FITTING (BAYESIAN MODEL) !             
!         +-----------------------------------------------+             
!                                                                       
!     THIS SUBROUTINE PRODUCES AN AUTOREGRESSIVE MODEL BY A BAYESIAN    
!     PROCEDURE USING THE OUTPUT OF SUBROUTINE REDUCT.                  
!                                                                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             MAICE                                                     
!             ARCOEF                                                    
!             BAYSPC                                                    
!             BAYSWT                                                    
!             COMAIC                                                    
!             PRINTA                                                    
!             SDCOMP                                                    
!       ----------------------------------------------------------------
!       INPUTS:                                                         
!         X:    (K+1)*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE REDU
!         D:    WORKING AREA                                            
!         K:    HIGHEST ORDER OF THE MODELS                             
!         N:    DATA LENGTH                                             
!         ISW:  =0     TO SUPPRESS THE OUTPUTS                          
!               =1     TO PRINT OUT THE OUTPUTS                         
!         TITLE:     TITLE OF DATA                                      
!         MJ1:  ABSOLUTE DIMENSION OF X                                 
!                                                                       
!       OUTPUTS:                                                        
!         A(I) (I=1,K):     AR-COEFFICIENTS                             
!         B(I) (I=1,K):     PARTIAL AUTOCORRELATION COEFFICIENTS        
!         SDB:  INNOVATION VARIANCE OF BAYESIAN MODEL                   
!         AICB: AIC OF BAYESIAN MODEL                                   
!                                                                       
!xx      IMPLICIT REAL * 8( A-H,O-Z )                                      
!c      REAL * 4   TITLE(20) , TTL(6)                                     
!c      DIMENSION  X(MJ1,1) , D(1) , A(1) , B(1)                          
!c      DIMENSION  SD(101) , AIC(101) , C(101)                            
!xx      DIMENSION  X(MJ1,1), D(K+1), A(K), B(K), B1(K)
!xx      DIMENSION  SD(K+1), AIC(K+1), C(K+1), DIC(K+1)
integer k, n, isw, mj1, imin
real(dp) x(mj1,1), d(k+1), sd(k+1), aic(k+1), dic(k+1),&
&aicm, sdmin, a(k), b1(k), b(k), c(k+1), sdb, pn,&
&aicb
! local
integer i
real(dp) fn
!c      DATA       TTL / 4H..  ,4H BAY,4HESIA,4HN MO,4HDEL ,4H  .. /      
!                                                                       
!c      IF(ISW.GE.1)  WRITE( 6,2 )                                        
fn = n
!                                                                       
!         +-----------------------------------------+                   
!         ! INNOVATION VARIANCE AND AIC COMPUTATION !                   
!         +-----------------------------------------+                   
!                                                                       
call  comaic( x,n,k,mj1,sd,aic )
!                                                                       
!         +-------------+                                               
!         ! AIC DISPLAY !                                               
!         +-------------+                                               
!                                                                       
!c      CALL  MAICE( AIC,SD,K,ISW,AICM,SDMIN,IMIN )                       
!xx      CALL MAICE( AIC,SD,K,ISW,AICM,SDMIN,IMIN,DIC )                   
call maice( aic,sd,k,aicm,sdmin,imin,dic )
!                                                                       
!         +-----------------------------+                               
!         ! BAYSIAN WEIGHTS COMPUTATION !                               
!         +-----------------------------+                               
!                                                                       
call  bayswt( aic,aicm,k,0,c )
!                                                                       
!         +-----------------------------------------------------+       
!         ! PARTIAL AUTOCORRELATION ESTIMATION (BAYESIAN MODEL) !       
!         +-----------------------------------------------------+       
!                                                                       
!c      CALL  BAYSPC( X,C,N,K,ISW,MJ1,B,D )                               
call  bayspc( x,c,n,k,isw,mj1,b1,b,d )
!                                                                       
!         +-----------------------------------------+                   
!         ! AUTOREGRESSIVE COEFFICIENTS COMPUTATION !                   
!         +-----------------------------------------+                   
!                                                                       
call  arcoef( b,k,a )
!                                                                       
!                                                                       
!          EQUIVALENT NUMBER OF PARAMETERS                              
!                                                                       
pn = 1.d0
do   10    i=1,k
!xx   10 PN = PN + D(I)*D(I)                                               
pn = pn + d(i)*d(i)
10 continue
!                                                                       
!         +-------------------------------------------+                 
!         ! INNOVATION VARIANCE OF THE BAYESIAN MODEL !                 
!         +-------------------------------------------+                 
!                                                                       
!c      CALL  SDCOMP( X,A,D,N,K,MJ1,SDB )                                 
call  sdcomp( x,a,n,k,mj1,sdb )
!                                                                       
!          -------------------------                                    
!          AIC OF THE BAYESIAN MODEL                                    
!          -------------------------                                    
!                                                                       
aicb = fn * dlog( sdb ) + 2.d0*pn
!                                                                       
!c      IF( ISW .EQ. 0 )     RETURN                                       
!                                                                       
!c      WRITE( 6,12 )                                                     
!c      WRITE( 6,9 )     (B(I),I=1,K)                                     
!c      WRITE( 6,15 )     SDB , PN , AICB                                 
!c      L1 = LAG + 1                                                      
!c      NL = N + LAG                                                      
!c      IF( ISW .GE. 1 )   CALL  PRINTA( A,SDB,K,TTL,6,TITLE,L1,NL )      
!                                                                       
return
!                                                                       
!cx    2 FORMAT( //1H ,18(1H-),/,' BAYESIAN PROCEDURE',/,1H ,18(1H-) )     
!cx    5 FORMAT( 1H ,4X,'I',10X,'SD',19X,'AIC',11X,'WEIGHTS' )             
!cx    6 FORMAT( 1H ,I5,D20.10,F16.3,F16.5 )                               
!cx    7 FORMAT( 1H ,'*****  MINIMUM AIC =',D17.10,3X,'ATTAINED AT M =',I3,
!cx     1'  *****' )                                                       
!cx    8 FORMAT( 1H ,'M =',I3,5X,'WEIGHT =',F8.5 )                         
!cx    9 FORMAT( 1H ,10F13.8 )                                             
!cx   11 FORMAT( //,1H ,132(1H-),/,' AR-COEFFICIENTS  (UP TO WEIGHT = 0.01)
!cx     1' )                                                               
!cx   12 FORMAT( ///' < BAYESIAN MODEL >',/,9H PARCOR'S )                  
!cx   13 FORMAT( ' AR-COEFFICIENTS' )                                      
!cx   15 FORMAT( 1H ,'INNOVATION VARIANCE',13X,'=',D24.8,/,1H ,'EQUIVALENT 
!cx     1NUMBER OF PARAMETERS =',F15.3,/,' EQUIVALENT AIC',18X,'=',F15.3 ) 
!cx   16 FORMAT( 1H+,32X,'***  MAICE MODEL  ***' )                         
!                                                                       
end
!
!
subroutine  arcoef( b,k,a )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PRODUCES THE AUTOREGRESSIVE COEFFICIENTS FROM A SE
!     OF PARTIAL AUTOCORRELATION COEFFICIENTS.                          
!                                                                       
!       INPUTS:                                                         
!         B:    VECTOR OF PARTIAL AUTOCORRELATIONS                      
!         K:    ORDER OF THE MODEL                                      
!                                                                       
!       OUTPUTS:                                                        
!         A:    VECTOR OF AR-COEFFICIENTS                               
!                                                                       
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  B(1) , A(1) , AA(100)                                  
!xx      DIMENSION  B(K) , A(K) , AA(K)
integer k
real(dp) b(k), a(k), aa(k)
! local
integer ii, im1, j, jj
!                                                                       
do  30     ii=1,k
a(ii) = b(ii)
aa(ii) = b(ii)
im1 = ii - 1
if( im1 .le. 0 )     go to 30
do  10     j=1,im1
jj = ii - j
!xx   10 A(J) = AA(J) - B(II)*AA(JJ)                                       
a(j) = aa(j) - b(ii)*aa(jj)
10 continue
if( ii .eq. k )     go to 40
do  20     j=1,im1
!xx   20 AA(J) = A(J)
aa(j) = a(j)
20 continue
30 continue
40 continue
return
end
!
!
!c      SUBROUTINE  ARMFIT( X,K,LAG,N,ISW,TITLE,MJ1,A,SDMIN,IMIN )        
subroutine  armfit( x,k,lag,n,isw,mj1,a,imin,sd,aic,dic,sdmin,&
!x     *                    AICM,IFG,LU )
&aicm )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          +------------------------------+                             
!          ! AUTOREGRESSIVE MODEL FITTING !                             
!          +------------------------------+                             
!                                                                       
!       THIS SUBROUTINE PRODUCES PARAMETERS OF AR-MODELS USING THE OUTPU
!       SUBROUTINE REDUCT.                                              
!                                                                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             MAICE                                                     
!             COMAIC                                                    
!             PRINTA                                                    
!             RECOEF                                                    
!       ----------------------------------------------------------------
!       INPUTS:                                                         
!         X:    (K+1)*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE REDU
!         K:    NUMBER OF REGRESSORS OF THE MODEL                       
!         LAG:  MAXIMUM TIME LAG OF THE MODEL                           
!         N:    DATA LENGTH                                             
!         ISW:  =0   TO PRODUCE THE MAICE MODEL ONLY (OUTPUTS SUPPRESSED
!               =1   TO PRODUCE THE MAICE MODEL ONLY                    
!               =2   TO PRODUCE ALL AR-MODELS (UP TO THE ORDER K)       
!         TITLE:     TITLE OF DATA                                      
!         MJ1:  ABSOLUTE DIMENSION OF X                                 
!                                                                       
!       OUTPUTS:                                                        
!         A:      VECTOR OF MAICE AR-COEFFICIENTS                       
!         SDMIN:  MAICE INNOVATION VARIANCE                             
!         IMIN:   MAICE ORDER                                           
!                                                                       
!xx      IMPLICIT REAL * 8( A-H,O-Z )                                      
!c      REAL * 4   TITLE(20) , TTL(4)                                     
!c      DIMENSION  X(MJ1,1) , A(1) , SD(101) , AIC(101)                   
!xx      DIMENSION  X(MJ1,1) , A(K)
!xx      DIMENSION  SD(K+1), AIC(K+1), DIC(K+1)
integer k, lag, n, isw, mj1, imin
real(dp) x(mj1,1), a(k), sd(k+1), aic(k+1), dic(k+1),&
&sdmin, aicm
! local
integer m, m1, l1, nl
!c      DATA     TTL / 4H   M,4H A I ,4H C E,4H   . /                    
!
!c      IF(ISW.GE.1)  WRITE( 6,3 )                                        
!x      IF( (ISW.GE.2) .AND. (IFG.NE.0) )  WRITE( LU,3 )
!                                                                       
!          +-----------------------------------------+                +-
!          ! INNOVATION VARIANCE AND AIC COMPUTATION !                ! 
!          +-----------------------------------------+                +-
!                                                                       
call  comaic( x,n,k,mj1,sd,aic )
!                                                                       
!          +-------------+                                            +-
!          ! AIC DISPLAY !                                            ! 
!          +-------------+                                            +-
!                                                                       
!c      CALL  MAICE( AIC,SD,K,ISW,AICM,SDMIN,IMIN )                       
!xx      CALL  MAICE( AIC,SD,K,ISW,AICM,SDMIN,IMIN,DIC )
call  maice( aic,sd,k,aicm,sdmin,imin,dic )
!                                                                       
!          +-----------------------------------------+                +-
!          ! AUTOREGRESSIVE COEFFICIENTS COMPUTATION !                ! 
!          +-----------------------------------------+                +-
!                                                                       
if( isw .lt. 2 )     go to  20
!                                                                       
!c           WRITE( 6,5 )                                                 
!x      IF ( IFG.NE.0 )  WRITE( LU,5 )
m = 0
!c      WRITE( 6,7 )     M , SD(1)                                        
!x      IF ( IFG.NE.0 )  WRITE( LU,7 )  M , SD(1)
do  10     m=1,k
m1 = m + 1
call  recoef( x,m,k,mj1,a )
!c   10 WRITE( 6,6 )     M , SD(M1) , (A(I),I=1,M)                        
!x   10 IF ( IFG.NE.0 )  WRITE( LU,6 )  M, SD(M1), (A(I),I=1,M)
10 continue
!                                                                       
20 if( imin .ge. 1 )  call  recoef( x,imin,k,mj1,a )
l1 = lag + 1
nl = n + lag
!c      IF( ISW .GE. 1 )   CALL  PRINTA( A,SDMIN,IMIN,TTL,4,TITLE,L1,NL ) 
!                                                                       
return
!                                                                       
!c    3 FORMAT( // 1H ,15(1H-),/,' MAICE PROCEDURE',/,1H ,15(1H-))        
!c    5 FORMAT( //' AR-COEFFICIENTS' )                                    
!c    6 FORMAT( ' M =',I3,9X,'(SD(M) =',D15.8,1X,')',/,(1X,10F13.8) )     
!c    7 FORMAT( ' M =',I3,9X,'(SD(M) =',D15.8,1X,')' )                    
!                                                                       
end
!
!c      SUBROUTINE  BAYSPC( X,C,N,K,ISW,MJ1,B,D )                         
subroutine  bayspc( x,c,n,k,isw,mj1,b,b1,d )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PRODUCES PARTIAL AUTOCORRELATION COEFFICIENTS B(I)
!     (I=1,K) OF THE BAYESIAN MODEL.                                    
!                                                                       
!         +-----------------------------------------------------+       
!         ! PARTIAL AUTOCORRELATION ESTIMATION (BAYESIAN MODEL) !       
!         +-----------------------------------------------------+       
!                                                                       
!       INPUTS:                                                         
!         X:    N*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE REDUCT  
!         C:    BAYESIAN WEIGHTS, OUTPUT OF SUBROUTINE BAYSWT           
!         N:    DATA LENGTH                                             
!         K:    HIGHEST ORDER OF THE MODELS                             
!         ISW:  =0     OUTPUTS ARE SUPPRESSED                           
!         ISW:  >0     OUTPUTS ARE PRINTED OUT                          
!         MJ1:  ABSOLUTE DIMENSION OF X                                 
!                                                                       
!       OUTPUTS:                                                        
!         B:    VECTOR OF PARTIAL AUTOCORRELATIONS                      
!         D:    INTEGRATED BAYESIAN WEIGHTS                             
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      REAL * 4  TABLE(41) , AST , AXIS , TERM                           
!x      CHARACTER TABLE(41),  AST , AXIS , TERM                           
!x      DIMENSION  X(MJ1,1) , B(1) , C(1) , D(1)                          
!x      DIMENSION  B1(1)
!xx      DIMENSION  X(MJ1,K+1) , B(K) , C(K+1) , D(K+1)                          
!xx      DIMENSION  B1(K)
integer n, k, isw, mj1
real(dp) x(mj1,k+1), c(k+1), b(k), b1(k), d(k+1)
! local
integer i, isig, isigm1, isigm2, isigp1, isigp2, j, j0, k1
real(dp) fn, sum, g, sc, sig
!c      DATA  AST , AXIS / 1H* , 1H! /                                    
!c      DATA  TABLE / 41*1H /                                             
!x      DATA  AST , AXIS / '*' , '!' /                                    
!x      DATA  TABLE/ 41*' ' /                                             
k1 = k + 1
fn = n
!                                                                       
!          PARTIAL AUTOCORRELATION (LEAST SQUARES ESTIMATE)             
!                                                                       
sum = x(k1,k1)**2
do 10  i=1,k
j = k1 - i
sum = sum + x(j,k1)**2
g = dsqrt( sum )
!xx   10 B(J) = X(J,K1)*X(J,J) / (G*DABS(X(J,J)))                          
b(j) = x(j,k1)*x(j,j) / (g*dabs(x(j,j)))
10 continue
!                                                                       
!          INTEGRATED BAYESIAN WEIGHT                                   
!                                                                       
d(k) = c(k1)
do  80     i=2,k
j = k1 - i
!xx   80 D(J) = D(J+1) + C(J+1)                                            
d(j) = d(j+1) + c(j+1)
80 continue
if( isw .eq. 0 )   go to 100
!                                                                       
!          PARCOR AND BAYESIAN WEIGHTS DISPLAY                          
!                                                                       
sc = 20.d0
sig = 1.d0/dsqrt(fn)
!xx      ISIG = SIG*SC + 0.5D0                                             
!xx      J0 = SC + 1                                                       
isig = int(sig*sc + 0.5d0)
j0 = int(sc) + 1
isigp1 = j0 + isig
isigm1 = j0 - isig
isigp2 = j0 + isig*2
isigm2 = j0 - isig*2
!x      TABLE(ISIGP1) = AXIS                                              
!x      TABLE(ISIGM1) = AXIS                                              
!x      TABLE(ISIGP2) = AXIS                                              
!x      TABLE(ISIGM2) = AXIS                                              
!x      TABLE(J0) = AXIS                                                  
!                                                                       
!c      WRITE( 6,11 )                                                     
!c      WRITE( 6,12 )                                                     
!x      DO  90     I=1,K                                                  
!x      BB = B(I)*SC                                                      
!x      IB = BB
!x      IF( BB .GT. 0.D0 )     IB = BB + 0.5D0                            
!x      IF( BB .LT. 0.D0 )     IB = BB - 0.5D0                            
!x      IB = IB + J0                                                      
!x      TERM = TABLE(IB)                                                  
!x      TABLE(IB) = AST                                                   
!c      WRITE( 6,13 )     I , B(I) , C(I+1) , D(I)                        
!c      WRITE( 6,14 )     (TABLE(J),J=1,41)                               
!x      TABLE(IB) = TERM                                                  
!x   90 CONTINUE                                                          
!                                                                       
!          PARCOR OF BAYESIAN MODEL                                     
!                                                                       
100 continue
do  110     i=1,k
!c  110 B(I) = B(I) * D(I)                                                
!xx  110 B1(I) = B(I) * D(I)
b1(i) = b(i) * d(i)
110 continue
!                                                                       
return
!                                                                       
!c   11 FORMAT( ///,73X,'PARCOR (LINES SHOW +SD AND +2SD), SD=SQRT(1/N)',/
!c     1,1H+,91X,'_',7X,'_' )                                             
!c   12 FORMAT( 1H ,12X,'PARCOR',9X,'WEIGHT',7X,'INTEGRATED',17X,'-1',19X,
!c     1'0',19X,'1',/,1H ,4X,'I',8X,'B(I)',11X,'W(I)',10X,'WEIGHT',19X,'+'
!c     2,4(10H---------+) )                                               
!c   13 FORMAT( 1H ,I5,3F15.8 )                                           
!c   14 FORMAT( 1H ,67X,41A1 )                                            
!cx   11 FORMAT(//65X,'PARCOR (LINES SHOW +/-SD AND +/-2SD), SD=SQRT(1/N)')
!cx   12 FORMAT( 1H ,12X,'PARCOR',9X,'WEIGHT',7X,'INTEGRATED',16X,'-1',19X,
!cx     1'0',19X,'1',/,1H ,4X,'I',8X,'B(I)',11X,'W(I)',10X,'WEIGHT',19X,'+'
!cx     2,4(10H---------+) )                                               
!cx   13 FORMAT( 1H ,I5,3F15.8,17X,41A1 )
!                                                                       
end
!
subroutine  bayswt( aic,aicm,k,isw,c )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE COMPUTES BAYESIAN WEIGHT OF AR-MODEL OF EACH ORDER
!                                                                       
!         +-----------------------------+                               
!         ! BAYESIAN WEIGHT COMPUTATION !                               
!         +-----------------------------+                               
!                                                                       
!     BAYESIAN WEIGHT OF THE M-TH ORDER MODEL IS DEFINED BY             
!            W(M)  =  CONST * C(M) / ( M+1 ),                           
!     WHERE                                                             
!            CONST =  NORMALIZING FACTOR                                
!            C(M)  =  EXP( -0.5*AIC(M) )                                
!                                                                       
!       INPUTS:                                                         
!         AIC:  VECTOR OF AIC'S                                         
!         AICM: MINIMUM AIC                                             
!         K:    HIGHEST ORDER OF THE MODELS                             
!         ISW:  =0   NON-ADAPTIVE DAMPER                                
!               =1   DATA ADAPTIVE DAMPER                               
!               =2   DAMPER DOES NOT USED                               
!                                                                       
!       OUTPUT:                                                         
!         C(I) (I=1,K+1):   VECTOR OF BAYESIAN WEIGHTS                  
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  AIC(1) , C(1)                                          
!xx      DIMENSION  AIC(K+1) , C(K+1)
integer k, isw
real(dp) aic(k+1), aicm, c(k+1)
! local
integer i, k1
real(dp) sum, ek, dic, ai
!                                                                       
!          C(I)  =  EXP( -AIC(I)/2 )                                    
!                                                                       
k1 = k + 1
sum = 0.d0
ek = 0.d0
!xx      DO  10     I=1,K1                                                 
do  11     i=1,k1
dic = -0.5d0 * (aic(i) - aicm)
c(i) = 0.d0
if( dic .lt. -40.d0 )     go to 10
c(i) = dexp(dic)
ek = ek + (i-1) * c(i)
10 sum = sum + c(i)
11 continue
!                                                                       
!          DAMPING OF C(I)                                              
!                                                                       
if( isw .eq. 1 )     go to 30
if( isw .eq. 2 )     go to 50
sum = 0.d0
do  20     i=1,k1
ai = i
c(i) = c(i) / ai
!xx   20 SUM = SUM + C(I)
sum = sum + c(i)
20 continue
go to 50
!                                                                       
!          DATA ADAPTIVE DAMPER                                         
!                                                                       
30 ek = ek / ( sum+ek )
sum = 0.d0
do  40     i=1,k1
c(i) = c(i) * ek**(i-1)
!xx   40 SUM = SUM + C(I)                                                  
sum = sum + c(i)
40 continue
!                                                                       
!          NORMALIZATION OF C                                           
!                                                                       
50 do  60     i=1,k1
!xx   60 C(I) = C(I) / SUM                                                 
c(i) = c(i) / sum
60 continue
!                                                                       
return
!                                                                       
end
!
!
subroutine  coef2( a,m,id,ii,jnd,lmax,mm,ksw,msw,mj1,b,c,e )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     COMPOSITION OF AR-COEFFICIENT MATRICES WITH INSTANTANEOUS RESPONSE
!                                                                       
!       INPUTS:                                                         
!          A(I)  (I=1,M):     REGRESSION COEFFICIENTS                   
!          JND(I) (I=1,M):    SPECIFICATION OF REGRESSORS               
!          M:     NUMBER OF REGRESSORS                                  
!          ID:    DIMENSION OF VECTOR OF OBSERVATIONS                   
!          II:    REGRESSAND                                            
!          LMAX:  HIGHEST ORDER OF THE MODEL                            
!          MM:    ORDER                                                 
!          KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR    
!                 =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR
!          MSW:   =0  CONSTANT VECTOR IS NOT ADOPTED AS A REGRESSOR     
!                 =1  CONSTANT VECTOR IS ACTUALLY ADOPTED AS REGRESSOR  
!          MJ1:   ABSOLUTE DIMENSION OF B                               
!       OUTPUTS:                                                        
!          B:     AR-COEFFICIENT MATRICES                               
!          C:     CONSTANT VECTOR                                       
!          E:     COEFFICIENT FOR INSTANTANEOUS RESPONSE                
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  A(1) , JND(1) , C(1)                                   
!x      DIMENSION  B(MJ1,MJ1,1) , E(MJ1,1)                                
!xx      DIMENSION  A(M) , JND(M) , C(ID)                                   
!xx      DIMENSION  B(MJ1,MJ1,MM) , E(MJ1,ID)                                
integer m, id, ii, jnd(m), lmax, mm, ksw, msw, mj1
real(dp) a(m), b(mj1,mj1,mm), c(id), e(mj1,id)
! local
integer i, j, jj, l, l1, m0
!                                                                       
m0 = msw + 1
c(ii) = 0.d0
if( msw .eq. 1 )     c(ii) = a(1)
do  100   jj=m0,m
i = jnd(jj) - ksw
l = i / id
j = i - l*id
if( j .ne. 0 )     go to 10
j = id
l = l - 1
10 l1 = l + 1
if( i .le. mm*id )     go to  20
e(ii,j) = -a(jj)
go to 30
20 b(ii,j,l1) = a(jj)
if( lmax .lt. l1 )     lmax = l1
30 continue
100 continue
do  40     i=1,id
!xx   40 E(I,I) = 1.D0                                                     
e(i,i) = 1.d0
40 continue
return
end
!
!
subroutine  comaic( x,n,k,mj1,sd,aic )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          +-----------------------------------------+                  
!          ! INNOVATION VARIANCE AND AIC COMPUTATION !                  
!          +-----------------------------------------+                  
!                                                                       
!       INPUTS:                                                         
!          X:     (K+1)*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE  R
!          N:     DATA LENGTH                                           
!          K:     NUMBER OF REGRESSORS                                  
!          MJ1:   ABSOLUTE DIMENSION OF X                               
!                                                                       
!       OUTPUTS:                                                        
!          SD(M)  (M=1,K+1):   VECTOR OF INNOVATION VARIANCES           
!          AIC(M) (M=1,K+1):   VECTOR OF AIC'S.                         
!                                                                       
!xx      IMPLICIT  REAL * 8( A-H,O-Z )                                     
!x      DIMENSION  X(MJ1,1) , AIC(1) , SD(1)                              
!xx      DIMENSION  X(MJ1,K+1) , AIC(K+1) , SD(K+1)                              
integer n, k, mj1
real(dp) x(mj1,k+1), sd(k+1), aic(k+1)
! local
integer i, k1, m
real(dp) fn, osd
fn = n
k1 = k + 1
!                                                                       
osd = 0.0d00
do 10   i = 1,k1
m = k1 - i + 1
osd = osd + x(m,k1)**2
sd(m)  = osd / fn
!xx   10 AIC(M) = FN*DLOG( SD(M) ) + 2*M                                   
aic(m) = fn*dlog( sd(m) ) + 2*m
10 continue
!                                                                       
return
!                                                                       
end
!
!
subroutine  copy( x,k,ii,jj,mj1,mj2,y )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!         +-----------------------+                                     
!         ! MAKE A COPY OF X ON Y !                                     
!         +-----------------------+                                     
!                                                                       
!        INPUTS:                                                        
!             X(I+II,J):     K*K MATRIX (I,J=1,...,K)                   
!             II:            INDICATES ORIGIN OF X                      
!             JJ:            INDICATES ORIGIN OF Y                      
!             MJ1:           ABSOLUTE DIMENSION OF X                    
!             MJ2:           ABSOLUTE DIMENSION OF Y                    
!                                                                       
!        OUTPUT:                                                        
!             Y(I+JJ,J):     COPY OF X (I,J=1,...,K)                    
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  X(MJ1,1) , Y(MJ2,1)                                    
!xx      DIMENSION  X(MJ1,K) , Y(MJ2,K)                                    
integer k, ii, jj, mj1, mj2
real(dp) x(mj1,k), y(mj2,k)
! local
integer i, i1, i2, j
!                                                                       
!xx      DO  10     I=1,K                                                  
do  20     i=1,k
i1 = i + ii
i2 = i + jj
do  10     j=1,k
!xx   10 Y(I2,J) = X(I1,J)
y(i2,j) = x(i1,j)
10 continue
20 continue
!                                                                       
return
!                                                                       
end
!
!
!c      SUBROUTINE  DELETE( X,D,IND,JND,K,L,M,MJ )                        
subroutine  delete( x,ind,jnd,k,l,m,mj )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!         +------------------------------------------------+            
!         ! DELETION OF THE VARIABLE M FROM THE REGRESSORS !            
!         +------------------------------------------------+            
!                                                                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
!             HUSHL1                                                    
!       ----------------------------------------------------------------
!       INPUTS:                                                         
!          IND(J)=I:     PRESENT STATUS,   VARIABLE J IS THE I-TH REGRES
!          JND(I)=J:     REQUIRED STATUS,  I-TH REGRESSOR IS VARIABLE J.
!          X:            (K+1)*(K+1)  MATRIX                            
!          D:            WORKING AREA                                   
!          K:            NUMBER OF VARIABLES                            
!          L:            NUMSER OF REGRESSORS IN THE PRESENT MODEL      
!          M:            INDICATION OF THE VARIABLE TO BE DELETED       
!          MJ:           ABSOLUTE DIMENSION OF X                        
!       OUTPUTS:                                                        
!          X:            (K+1)*(K+1)  MATRIX                            
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  X(MJ,1) , D(1) , IND(1) , JND(1)                       
!x      DIMENSION  X(MJ,1) , IND(1) , JND(1)                       
!xx      DIMENSION  X(MJ,K+1) , IND(K+1) , JND(K+1)
integer k, l, m, mj, ind(k+1), jnd(k+1)
real(dp) x(mj,k+1)
! local
integer i, ii, i1, j, k1, lm1
!
k1 = k + 1
do  60     i=1,k1
j = jnd(i)
!xx   60 IND(J) = I                                                        
ind(j) = i
60 continue
ii = ind(m)
!c      IF( II-L )     10,30,40                                           
if( ii-l .lt. 0 )  go to 10
if( ii-l .eq. 0 )  go to 30
if( ii-l .gt. 0 )  go to 40
10 i1 = ii + 1
do  20     i=i1,l
!xx   20 JND(I-1) = JND(I)                                                 
jnd(i-1) = jnd(i)
20 continue
jnd(l) = m
lm1 = l - 1
!                                                                       
!c      CALL  HUSHL1( X,D,MJ,K1,LM1,II,IND,JND )                          
call  hushl1( x,mj,k1,lm1,ii,ind,jnd )
!                                                                       
30 l = l - 1
!                                                                       
40 return
end
!
!
subroutine fouger(g,lgp1,fc,fs,lf1)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     FOURIER TRANSFORM (GOERTZEL METHOD)                               
!     THIS SUBROUTINE COMPUTES FOURIER TRANSFORM OF G(I),I=0,1,...,LG AT
!     FREQUENCIES K/(2*LF), K=0,1,...,LF AND RETURNS COSINE TRANSFORM IN
!     FC(K) AND SINE TRANSFORM IN FS(K).                                
!                                                                       
!       INPUTS:                                                         
!          G(I):   ORIGINAL DATA (I=0,1,...,LG)                         
!          LG1:    = LG+1                                               
!          LF1:    = LF+1                                               
!                                                                       
!       OUTPUTS:                                                        
!          FC(I):  COSINE TRANSFORM OF G  (I=0,1,...,LF)                
!          FB(I):  SINE TRANSFORM OF G  (I=0,1,...,LF)                  
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!xx      DIMENSION G(LGP1),FC(LF1),FS(LF1)
integer lgp1, lf1
real(dp) g(lgp1), fc(lf1), fs(lf1)
! local
integer i, i2, k, lg, lf, lg3, lg4
real(dp) t, pi, alf, ak, tk, ck, sk, ck2, um0, um1, um2
lg=lgp1-1
lf=lf1-1
!     REVERSAL OF G(I),I=1,...,LGP1 INTO G(LG3-I)   LG3=LGP1+1          
if(lgp1.le.1) go to 110
lg3=lgp1+1
lg4=lgp1/2
do 100 i=1,lg4
i2=lg3-i
t=g(i)
g(i)=g(i2)
!xx  100 G(I2)=T                                                           
g(i2)=t
100 continue
110 pi=3.1415926536d0
alf=lf
t=pi/alf
do 10 k=1,lf1
ak=k-1
tk=t*ak
ck=dcos(tk)
sk=dsin(tk)
ck2=ck+ck
um1=0.0d0
um2=0.0d0
if(lg.eq.0) go to 12
do 11 i=1,lg
um0=ck2*um1-um2+g(i)
um2=um1
!xx   11 UM1=UM0
um1=um0
11 continue
12 fc(k)=ck*um1-um2+g(lgp1)
fs(k)=sk*um1
10 continue
return
end
!
!
!c      SUBROUTINE  HUSHL1( X,D,MJ1,K,L,M,IND,JND )                       
subroutine  hushl1( x,mj1,k,l,m,ind,jnd )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PERFORMS THE HOUSEHOLDER TRANSFORMATION OF THE MAT
!                                                                       
!       INPUTS:                                                         
!         X:   ORIGINAL (K+1)*(K+1) MATRIX                              
!         N:   NUMBER OF ROWS OF X,  NOT GREATER THAN MJ1               
!         K:   NUMBER OF COLUMNS OF X                                   
!         L:   END POSITION OF THE HOUSEHOLDER TRANSFORMATION           
!         M:   STARTING POSITION OF THE HOUSEHOLDER TRANSFORMATION      
!         IND:   SPECIFICATION OF THE PRESENT FORM OF X                 
!                IND(J) = I;     VARIABLE I IS THE J-TH REGRESSOR       
!         JND:   SPECIFICATION OF THE REQUIRED FORM OF X                
!                JND(I) = J;     THE I-TH REGRESSOR IS VARIABLE J       
!                                                                       
!       OUTPUTS:                                                        
!         X:   TRANSFORMED MATRIX                                       
!                                                                       
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  X(MJ1,1) , D(MJ1) , IND(1) , JND(1)                    
!xx      DIMENSION  X(MJ1,K) , D(MJ1) , IND(K) , JND(K)
integer mj1, k, l, m, ind(k), jnd(k)
real(dp) x(mj1,k)
! local
integer i, ii, ii1, j, jj, j1, nn
real(dp) d(mj1), tol, h, g, f, s
!                                                                       
tol = 1.0d-60

!
nn = 0
!xx      DO  100     II=M,L                                                
do  110     ii=m,l
jj = jnd(ii)
nn = max0( nn,ind(jj) )
h = 0.0d00
do  10     i=ii,nn
d(i) = x(i,jj)
!xx   10 H = H + D(I)*D(I)                                                 
h = h + d(i)*d(i)
10 continue
if( h .gt. tol )     go to 20
g = 0.0d00
go to 100
20 g = dsqrt( h )
f = x(ii,jj)
if( f .ge. 0.0d00 )     g = -g
d(ii) = f - g
h = h - f * g
!                                                                       
!     ( I - D*D'/H ) * X                                                
!                                                                       
ii1 = ii + 1
if( ii .eq. nn )   go to 35
do  30     i=ii1,nn
!xx   30 X(I,JJ) = 0.D0                                                    
x(i,jj) = 0.d0
30 continue
35 continue
if( ii .eq. k )     go to 100
do  60     j1=ii1,k
j = jnd(j1)
s = 0.0d00
do  40     i=ii,nn
!xx   40 S = S + D(I)*X(I,J)
s = s + d(i)*x(i,j)
40 continue
s = s / h
do  50     i=ii,nn
!xx   50 X(I,J) = X(I,J) - D(I)*S
x(i,j) = x(i,j) - d(i)*s
50 continue
60 continue
100 x(ii,jj) = g
110 continue
return
end
!
!
!c      SUBROUTINE  HUSHLD( X,D,MJ1,N,K )                                 
subroutine  hushld( x,mj1,n,k )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          +----------------------------+                               
!          ! HOUSEHOLDER TRANSFORMATION !                               
!          +----------------------------+                               
!                                                                       
!     THIS SUBROUTINE  TRANSFORMS MATRIX X INTO AN UPPER TRIANGULAR FORM
!     BY HOUSEHOLDER TRANSFORMATION.                                    
!                                                                       
!       INPUTS:                                                         
!          X:     ORIGINAL N*K DATA MATRIX                              
!          D:     WORKING AREA                                          
!          MJ1:   ABSOLUTE DIMENSION OF X                               
!          N:     NUMBER OF ROWS OF X, NOT GREATER THAN MJ1             
!          K:     NUMBER OF COLUMNS OF X                                
!                                                                       
!       OUTPUT:                                                         
!          X:     SQUARE ROOT OF DATA COVARIANCE MATRIX (UPPER TRIANGULA
!                                                                       
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  X(MJ1,1) , D(MJ1)                                      
!xx      DIMENSION  X(MJ1,K) , D(MJ1)
integer mj1, n, k
real(dp) x(mj1,k)
! local
integer i, ii, ii1, j
real(dp) d(mj1), tol, h, g, f, s
!                                                                       
tol = 1.0d-60
!                                                                       
!xx      DO 100  II=1,K                                                    
do 110  ii=1,k
h = 0.0d00
do 10  i=ii,n
d(i) = x(i,ii)
!xx   10       H = H + D(I)*D(I)                                           
h = h + d(i)*d(i)
10 continue
if( h .gt. tol )  go to 20
g = 0.0d00
go to 100
20 g = dsqrt( h )
f = x(ii,ii)
if( f .ge. 0.0d00 )   g = -g
d(ii) = f - g
h = h - f*g
!                                                                       
!          FORM  (I - D*D'/H) * X, WHERE H = D'D/2                      
!                                                                       
ii1 = ii+1
do 30  i=ii1,n
!xx   30    X(I,II) = 0.0D00                                               
x(i,ii) = 0.0d00
30 continue
if( ii .eq. k )  go to 100
do 60  j=ii1,k
s = 0.0d00
do 40  i=ii,n
!xx   40       S = S + D(I)*X(I,J)                                         
s = s + d(i)*x(i,j)
40 continue
s = s/h
do 50  i=ii,n
!xx   50      X(I,J) = X(I,J) - D(I)*S                                     
x(i,j) = x(i,j) - d(i)*s
50 continue
60 continue
100 x(ii,ii) = g
110 continue
!                                                                       
return
!                                                                       
end
!

!
!c      SUBROUTINE  MAICE( AIC,SD,K,ISW,AICM,SDM,IMIN )                   
!xx      SUBROUTINE  MAICE( AIC,SD,K,ISW,AICM,SDM,IMIN,DIC )
subroutine  maice( aic,sd,k,aicm,sdm,imin,dic )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!             +-------------+                                           
!             ! AIC DISPLAY !                                           
!             +-------------+                                           
!                                                                       
!        THIS SUBROUTINE PRODUCES NUMERICAL AND GRAPHICAL DISPLAYS OF AI
!                                                                       
!       INPUTS:                                                         
!          AIC:   VECTOR OF AIC'S                                       
!          SD:    VECTOR OF INNOVATION VARIANCES                        
!          K:     UPPER LIMIT OF THE ORDER                              
!          ISW:   =0   OUTPUTS ARE SUPPRESSED                           
!                 >0   AIC'S ARE DIPLAIED                               
!                                                                       
!       OUTPUTS:                                                        
!          AICM:  MINIMUM AIC                                           
!          SDM:   MAICE INNOVATION VARIANCE                             
!          IMIN:  MAICE ORDER                                           
!                                                                       
!xx      IMPLICIT  REAL * 8( A-H,O-Z )                                     
!c      REAL * 4  TBL(41) , AST , BLNK , AXIS , PERI                      
!x      DIMENSION  AIC(1) , SD(1)                                         
!xx      DIMENSION  AIC(K+1) , SD(K+1)                                         
!xx      DIMENSION  DIC(K+1)
integer k, imin
real(dp) aic(k+1), sd(k+1), aicm, sdm, dic(k+1)
! local
integer i, k1
!c      DATA  AST / 1H* / , BLNK / 1H  / , AXIS / 1H! / , PERI / 1H. /    
!                                                                       
!c      DO  10     I=1,41                                                 
!c   10 TBL(I) = BLNK                                                     
!                                                                       
!       SEARCH FOR THE MINIMUM OF AIC(I)                                
!                                                                       
k1 = k + 1
imin = 0
sdm  = sd(1)
aicm = aic(1)
do  20   i = 1,k
if( aic(i+1) .ge. aicm )  go to 20
imin = i
sdm = sd(i+1)
aicm = aic(i+1)
20 continue
do 25    i = 1,k1
!xx   25 DIC(I) = AIC(I) - AICM
dic(i) = aic(i) - aicm
25 continue
!c      IF( ISW .LE. 0 )     RETURN                                       
!                                                                       
!       DISPLAY OF AIC'S                                                
!                                                                       
!c      WRITE( 6,5 )                                                      
!c      DO  30   I = 1,K1                                                 
!c      II = I - 1                                                        
!c      DIC = AIC(I) - AICM                                               
!c      TBL(1) = AXIS                                                     
!c      ID = DIC + 1.5D0                                                  
!c      ID = DIC(I) + 1.5D0                                                  
!c      IF( ID .LE. 41 )     TBL(ID) = AST                                
!c      IF( ID .GT. 41 )     TBL(41) = PERI                               
!c      WRITE( 6,6 )     II , SD(I) , AIC(I) , DIC , (TBL(J),J=1,41)      
!c      IF( ID .LE. 41 )     TBL(ID) = BLNK                               
!c      IF( ID .GT. 41 )     TBL(41) = BLNK                               
!c   30 CONTINUE                                                          
!c      WRITE( 6,7 )     AICM , IMIN , SDM                                
!                                                                       
return
!                                                                       
!cx    5 FORMAT( /1H ,70X,'AIC(M)-AICMIN (TRUNCATED AT 40.0)',/,2X,'ORDER',
!cx     1 2X,'INNOVATION VARIANCE',40X,'0',8X,'10',8X,'20',8X,'30',8X,'40',
!cx     2/,5X,'M',10X,'SD(M)',14X,'AIC(M)',7X,'AIC(M)-AICMIN',7X,'+',4(10H-
!cx     3--------+) )                                                      
!cx    6 FORMAT( 1H ,I5,D20.10,2F16.3,10X,41A1 )                           
!cx    7 FORMAT( /1H ,'*****  MINIMUM AIC =',D17.10,3X,'ATTAINED AT M =',
!cx     1I3,5X,'SD(M) =',D17.10,'  *****' )
!                                                                       
end
!
!
subroutine  marcof( d,e,id,m,mj3,a,b )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE COMPUTES COEFFICIENT MATRICES OF MULTI-VARIATE AUT
!     REGRESSIVE MODEL FROM PARTIAL AUTOREGRESSION COEFFICIENT MATRICES 
!     FORWARD AND BACKWARD MODELS.                                      
!                                                                       
!       INPUTS:                                                         
!          D:     PARTIAL AUTOCORRELATIONS OF FORWARD MODEL             
!          E:     PARTIAL AUTOCORRELATIONS OF BACKWARD MODEL            
!          ID:    DIMENSION OF THE PROCESS                              
!          M:     ORDER OF THE MODELS                                   
!          MJ3:   ABSOLUTE DIMENSION OF A,B,D AND E IN THE MAIN PROGRAM 
!                                                                       
!       OUTPUTS:                                                        
!          A:     AR-COEFFICIENT MATRICES OF FORWARD MODEL              
!          B:     AR-COEFFICIENT MATRICES OF BACKWARD MODEL             
!xx      IMPLICIT  REAL * 8 ( A-H , O-Z )                                  
!x      DIMENSION  A(MJ3,MJ3,1) , B(MJ3,MJ3,1)                            
!x      DIMENSION  D(MJ3,MJ3,1) , E(MJ3,MJ3,1)                            
!c      DIMENSION  F(10,10) , G(10,10)                                    
!xx      DIMENSION  A(MJ3,MJ3,M) , B(MJ3,MJ3,M)                            
!xx      DIMENSION  D(MJ3,MJ3,M) , E(MJ3,MJ3,M)                            
!xx      DIMENSION  F(ID,ID) , G(ID,ID)                                    
integer id, m, mj3
real(dp) d(mj3,mj3,m), e(mj3,mj3,m), a(mj3,mj3,m),&
&b(mj3,mj3,m)
! local
integer i, ii, ii1, imj, j, jj, l
real(dp) f(id,id), g(id,id), suma, sumb
!                                                                       
!xx      DO 10  II=1,M                                                     
!xx      DO 10  I=1,ID                                                     
do 12  ii=1,m
do 11  i=1,id
do 10  j=1,id
a(i,j,ii) = d(i,j,ii)
!xx   10 B(I,J,II) = E(I,J,II)                                             
b(i,j,ii) = e(i,j,ii)
10 continue
11 continue
12 continue
if(  m .eq. 1 )  return
!                                                                       
do 60  ii=2,m
!                                                                       
ii1 = ii - 1
do 50  jj=1,ii1
imj = ii - jj
!xx      DO 20  I=1,ID
do 21  i=1,id
do 20  j=1,id
f(i,j) = a(i,j,imj)
!xx   20 G(I,J) = B(I,J,JJ)
g(i,j) = b(i,j,jj)
20 continue
21 continue
!xx      DO 40  I=1,ID                                                     
do 41  i=1,id
do 40  j=1,id
suma = f(i,j)
sumb = g(i,j)
do 30  l=1,id
suma = suma - a(i,l,ii)*g(l,j)
!xx   30 SUMB = SUMB - B(I,L,II)*F(L,J)                                    
sumb = sumb - b(i,l,ii)*f(l,j)
30 continue
a(i,j,imj) = suma
!xx   40 B(I,J,JJ) = SUMB                                                  
b(i,j,jj) = sumb
40 continue
41 continue
50 continue
!                                                                       
60 continue
return
!                                                                       
end
!
!
!c      SUBROUTINE  MARFIT( X,Y,D,N,ID,M,KSW,MJ1,MJ2,MJ3,MJ4,ISW,IPR,B,E, 
!c     *                    EX,C,LMAX,AICS )                              
subroutine  marfit( x,n,id,m,ksw,mj1,mj2,mj3,mj4,isw,ipr,aic,sd,&
!x     *DIC,AICM,SDM,IM,BI,EI,B,E,EX,C,LMAX,AICS,JNDF,AF,NPR,AAIC,IFG,LU )
&dic,aicm,sdm,im,bi,ei,b,e,ex,c,lmax,aics,jndf,af,npr,aaic )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!         MULTI-VARIATE AUTOREGRESSIVE MODEL FITTING                    
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             COPY                                                      
!             COEF2                                                     
!             MAICE                                                     
!             MCOEF                                                     
!             ADDVAR                                                    
!             AICCOM                                                    
!             DELETE                                                    
!             HUSHL1                                                    
!             SRCOEF                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          X:     ((M+1)*ID)*((M+1)*ID)  UPPER TRIANGULAR MATRIX,       
!                 OUTPUT OF SUBROUTINE MREDCT                           
!          Y:     WORKING AREA (MATRIX)                                 
!          D:     WORKING AREA                                          
!          N:     DATA LENGTH                                           
!          ID:    DIMENSION OF DATA                                     
!          M:     HIGHEST ORDER OF THE MODELS                           
!          KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR    
!                 =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR
!          MJ1:   ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM           
!          MJ2:   ABSOLUTE DIMENSION OF E IN THE MAIN PROGRAM           
!          MJ3:   ABSOLUTE DIMENSION OF B IN THE MAIN PROGRAM           
!          MJ4:   ABSOLUTE DIMENSION OF Y IN THE MAIN PROGRAM           
!          ISW:   =0  MULTI-VARIATE AUTOREGRESSIVE MODEL IS REQUESTED   
!                 =1  INSTANTANEOUS RESPONSE MODEL IS REQUESTED         
!          IPR:   PRINT OUT CONTROL                                     
!                                                                       
!       OUTPUTS:                                                        
!          B:     AR-COEFFICIENT MATRICES                               
!          E:     INNOVATION VARIANCE MATRIX                            
!          EX:    RESIDUAL VARIANCES OF INSTANTENEOUS RESPONSE MODELS   
!          C:     CONSTANT VECTOR                                       
!          LMAX:  ORDER OF THE MAICE MODEL                              
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  X(MJ1,1) , Y(MJ4,1) , D(1) , E(MJ2,1) , B(MJ2,MJ2,MJ3) 
!c      DIMENSION  C(1) , A(100) , AIC(51) , SD(51) , EX(1)               
!c      DIMENSION  IND(100) , JND(100) , KND(100)                         
!x      DIMENSION  X(MJ1,1), Y(MJ4,MJ4), E(ID,ID), B(ID,ID,M) 
!xx      DIMENSION  X(MJ1,MJ4), Y(MJ4,MJ4), E(ID,ID), B(ID,ID,M) 
!xx      DIMENSION  C(ID), A(MJ4)
!xx      DIMENSION  AIC(M+1,ID), SD(M+1,ID), DIC(M+1,ID), EX(ID)
!xx      DIMENSION  AICM(ID), SDM(ID), IM(ID)
!xx      DIMENSION  IND(MJ4) , JND(MJ4) , KND(MJ4)                         
!xx      DIMENSION  EI(ID,ID) , BI(ID,ID,M)
!xx      DIMENSION  JNDF(MJ4,ID), AF(MJ4,ID)
!xx      DIMENSION  NPR(ID), AAIC(ID)
integer n, id, m, ksw, mj1, mj2, mj3, mj4, isw, ipr, im(id),&
&lmax, jndf(mj4,id), npr(id)
real(dp) x(mj1,mj4), aic(m+1,id), sd(m+1,id), dic(m+1,id),&
&aicm(id), sdm(id), bi(id,id,m), ei(id,id),&
&b(id,id,m), e(id,id), ex(id), c(id), aics,&
&af(mj4,id), aaic(id)
! local
integer i, i0, i1, ii, imin, imp1, ipr2, j, j0, j1, jj, jm1, k,&
&k0, k01, k02, kk, kk1, kx1, m1, madd, md, mdel, md0, md2,&
&msw, nsw, ind(mj4), jnd(mj4), knd(mj4)
real(dp) y(mj4,mj4), a(mj4), aicsum, osd, oaic, aicmin,&
&sdmin, sdd, aic1, aic2, osp
!                                     
!                                                                       
!         INITIAL SETTING                                               
!                                                                       
m1 = m + 1
md0 = m * id
md = m*id + ksw
md2 = m1*id + ksw
aicsum = 0.d0
lmax = 0
nsw = 0
!xx      DO 20  I=1,MJ2                                                    
!xx      DO 20  J=1,MJ2                                                    
!xx      DO 10  K=1,MJ3                                                    
!xx   10 B(I,J,K) = 0.D0                                                   
!xx   20 E(I,J) = 0.D0
b(1:mj2,1:mj2,1:mj3) = 0.d0
e(1:mj2,1:mj2) = 0.d0
do 30  i=1,md2
ind(i) = i
!xx   30 JND(I) = I                                                        
jnd(i) = i
30 continue
!                                                                       
call  copy( x,md2,0,md2,mj1,mj1,x )
!                                                                       
!                                                                       
do 500     ii=1,id
msw = ksw
!                                                                       
!c      IF(IPR.GE.3)  WRITE( 6,3 )                                        
!c      IF( IPR.GE.2)   WRITE( 6,645 )   II                               
!c      IF(IPR.GE.3)  WRITE( 6,642 )                                      
!x      IF( (IPR.GE.2) .AND. (IFG.NE.0) )  WRITE( LU,645 )   II
!x      IF( (IPR.GE.3) .AND. (IFG.NE.0) )  WRITE( LU,642 )
jj = ii - 1
kk = md + jj
kk1 = kk + 1
!                                                                       
!         ADDITION OF REGRESSOR (INSTANTANEOUS RESPONSE FROM II-TH VARIA
!                                                                       
call  copy( x,md2,md2,0,mj1,mj4,y )
call  copy( x,md2,md2,0,mj1,mj1,x )
!                                                                       
if( ii .le. 1 )  go to 40
j1 = jj + ksw
do 70  i=1,md2
j = jnd(i)
!xx   70 IND(J) = I                                                        
ind(j) = i
70 continue
jnd(1) = 1
do 75  i=1,jj
j = md + i
i1 = i + ksw
!xx   75 JND(I1) = J                                                       
jnd(i1) = j
75 continue
do 80  i=1,md0
j = i + ksw
i1 = i + ksw + jj
!xx   80 JND(I1) = J                                                       
jnd(i1) = j
80 continue
i1 = jj + 1
do 85  i=i1,id
j = md + i
!xx   85 JND(J) = J                                                        
jnd(j) = j
85 continue
!                                                                       
!c      CALL  HUSHL1( X,D,MJ1,MD2,MD2,1,IND,JND )                         
call  hushl1( x,mj1,md2,md2,1,ind,jnd )
!                                                                       
call  copy( x,md2,0,md2,mj1,mj1,x )
!                                                                       
!--------------------   FIRST STEP OF AIC MINIMIZATION   ---------------
!                                                                       
!          AIC'S OF INITIAL MODELS COMPUTATION                          
!                                                                       
40 do 50  i=1,m1
k = (i-1)*id + jj + ksw
call  aiccom( x,n,k,kk,mj1,osd,oaic )
!c      SD(I) = OSD                                                       
!c   50 AIC(I) = OAIC                                                     
sd(i,ii) = osd
!xx   50 AIC(I,II) = OAIC                                                  
aic(i,ii) = oaic
50 continue
do 60  i=1,md2
!xx   60 KND(I) = JND(I)                                                   
knd(i) = jnd(i)
60 continue
!                                                                       
!         ORDER DETERMINATION BY AIC ( INITIAL ESTMATE )                
!                                                                       
ipr2 = ipr - 2
!c      CALL  MAICE( AIC,SD,M,IPR2,AICMIN,SDMIN,IMIN )                    
!xx      CALL MAICE( AIC(1,II),SD(1,II),M,IPR2,AICMIN,SDMIN,IMIN,
!xx     *            DIC(1,II) )
call maice( aic(1,ii),sd(1,ii),m,aicmin,sdmin,imin,dic(1,ii) )
aicm(ii) = aicmin
sdm(ii) = sdmin
im(ii) = imin
!                                                                       
k0 = imin*id + jj + ksw
!                                                                       
!                                                                       
!         REGRESSION COEFFICIENTS COMPUTATION ( INITIAL ESTIMATE )      
!                                                                       
if( ipr .lt. 3 )  go to 90
!c      WRITE( 6,4 )                                                      
!c      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD )                       
!c      WRITE( 6,643 )                                                    
!x      IF ( IFG.NE.0 )  WRITE( LU,4 )
!x      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD,AAIC(II),IFG,LU )
if( imin.gt.0 )&
!xx     *    CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD,AAIC(II) )
&call  srcoef( x,k0,kk,n,mj1,jnd,a,sdd,aaic(ii) )
!x      IF ( IFG.NE.0 )  WRITE( LU,643 )
90 continue
!                                                                       
!                                                                       
!--------------------   SECOND STEP OF AIC MINIMIZATION   --------------
!                                                                       
imp1 = imin + 1
do 230  i1=1,id
do 100  i=1,kk
!xx  100 KND(I) = JND(I)                                                   
knd(i) = jnd(i)
100 continue
!                                                                       
call  copy( x,kk1,0,0,mj1,mj4,y )
!c      IF(IPR.GE.3)  WRITE( 6,11 )     I1                                
!x      IF( (IPR.GE.3) .AND. (IFG.NE.0) )  WRITE( LU,11 )     I1
aic1 = aicmin
aic2 = aicmin
k01 = k0
!c      IMP1 = IMIN + 1                                                   
!c      IF(IPR.GE.3)  WRITE( 6,13 )     AIC1                              
!x      IF( (IPR.GE.3) .AND. (IFG.NE.0) )  WRITE( LU,13 )     AIC1
if( imin .ge. m )  go to 140
!                                                                       
!          CHECK REGRESSOR MADD  < ADD? >                               
!                                                                       
do 110  j1=imp1,m
madd = (j1-1)*id + i1 + ksw
if( ind(madd) .eq. 1 )  go to 110
k = k01 + 1
kx1 = k
!c      CALL  ADDVAR( X,D,IND,JND,KK,KX1,MADD,MJ1 )                       
call  addvar( x,ind,jnd,kk,kx1,madd,mj1 )
call  aiccom( x,n,k,kk,mj1,osd,oaic )
!                                                                       
!         DECISION BY AIC                                               
!                                                                       
!c      IF(IPR.GE.3)  WRITE( 6,5 )   MADD,OAIC,MADD,J1,I1                 
if( oaic .gt. aic1 )  go to 120
!c      IF( IPR .GE. 3 )     WRITE( 6,6 )                                 
!x      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,56 ) MADD,OAIC,MADD,J1,I1
aic1 = oaic
k01 = k01 + 1
110 continue
go to 140
120 k = k - 1
!c      IF( IPR .GE. 3 )     WRITE( 6,9 )                                 
!x      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,59 ) MADD,OAIC,MADD,J1,I1
if( j1 .ne. imp1 )  go to 140
call  copy( y,kk1,0,0,mj4,mj1,x )
do 130  j=1,kk
i = ind(j)
!xx  130 JND(I) = J                                                        
jnd(i) = j
130 continue
140 continue
!                                                                       
!         CHECK REGRESSOR MDEL  < DELETE ? >                            
!                                                                       
if( imin .lt. 1 )  go to 200
k02 = k0
do 170  jm1=1,imin
j1 = imin - jm1 + 1
mdel = (j1-1)*id + i1 + ksw
k = k02 - 1
do 150  i0=1,md2
if( jnd(i0) .eq. mdel )  go to 160
150 continue
go to 170
160 continue
!c      CALL  DELETE( Y,D,IND,KND,KK,K02,MDEL,MJ4 )                       
call  delete( y,ind,knd,kk,k02,mdel,mj4 )
call  aiccom( y,n,k,kk,mj4,osd,oaic )
!                                                                       
!         DECISION BY AIC                                               
!                                                                       
!c      IF(IPR.GE.3)  WRITE( 6,7 )   MDEL,OAIC,MDEL,J1,I1                 
if( oaic .gt. aic2 )  go to 180
!c      IF( IPR .GE. 3 )     WRITE( 6,8 )                                 
!x      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,78 ) MDEL,OAIC,MDEL,J1,I1
aic2 = oaic
if( aic2 .ge. aic1 )  go to 170
do 165  i=1,kk
!xx  165 JND(I) = KND(I)                                                   
jnd(i) = knd(i)
165 continue
call  copy( y,kk1,0,0,mj4,mj1,x )
170 continue
go to  200
180 k02 = k02 + 1
!c      IF( IPR .GE. 3 )     WRITE( 6,18 )                                
!x      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,718 ) MDEL,OAIC,MDEL,J1,I1
200 continue
aicmin = dmin1( aic1,aic2 )
!                                                                       
!         COMPARISON OF AIC1 AND AIC2                                   
!                                                                       
k0 = k01
if( aic1 .le. aic2 )  go to 220
k0 = k02
220 continue
!c      IF(IPR.GE.3)  WRITE( 6,12 )     AICMIN ,K0                        
!x      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,12 )  AICMIN ,K0
230 continue
!                                                                       
!                                                                       
!       SECOND ESTIMATE                                                 
!                                                                       
if( ipr .lt. 3 )  go to 240
!c      WRITE( 6,14 )                                                     
!c      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD )                       
!c      WRITE( 6,644 )                                                    
!x      IF( IFG.NE.0 )  WRITE( LU,14 )
!x      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD,AAIC(II),IFG,LU )
!xx      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD,AAIC(II) )
call  srcoef( x,k0,kk,n,mj1,jnd,a,sdd,aaic(ii) )
!x      IF( IFG.NE.0 )  WRITE( LU,644 )                                   
240 continue
!                                                                       
!--------------------   FINAL STEP OF AIC MINIMIZATION   ---------------
!                                                                       
if(ksw.ne.1)  go to 280
!                                                                       
!          CHECK CONSTANT VECTOR  < DELETE ? >                          
!                                                                       
call  copy( x,kk1,0,0,mj1,mj4,y )
!c      IF(IPR.GE.3)  WRITE( 6,19 )                                       
!c      IF(IPR.GE.3)  WRITE( 6,13 )  AICMIN                               
!x      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,19 )
!x      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,13 )  AICMIN
mdel = 1
k = k0-1
!c      CALL  DELETE( Y,D,IND,JND,KK,K0,MDEL,MJ4 )                        
call  delete( y,ind,jnd,kk,k0,mdel,mj4 )
!                                                                       
call  aiccom( y,n,k,kk,mj4,osp,oaic )
!c      IF(IPR.GE.3)  WRITE( 6,21)  MDEL,OAIC,MDEL                        
if(oaic.ge.aicmin)  go to 250
!                                                                       
!c      IF(IPR.GE.3)  WRITE( 6,8 )                                        
!x      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,218 )  MDEL,OAIC,MDEL
aicmin = oaic
msw = 0
call  copy( y,kk1,0,0,mj4,mj1,x )
go to 270
!                                                                       
250 k0 = k0+1
!c      IF(IPR.GE.3)  WRITE( 6,18 )                                       
!x      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,2118 )  MDEL,OAIC,MDEL
do 260  i=1,kk
j = ind(i)
!xx  260 JND(J) = I                                                        
jnd(j) = i
260 continue
!c  270 IF(IPR.GE.3)  WRITE( 6,16 )  AICMIN,K0                            
!x  270 IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,16 )  AICMIN,K0
270 continue
280 continue
!                                                                       
!         CHECK REGRESSOR MDEL  < DELETE ? >                            
!                                                                       
do 400     i1=1,id
!                                                                       
call  copy( x,kk1,0,0,mj1,mj4,y )
!                                                                       
!c      IF(IPR.GE.3)  WRITE( 6,11 )     I1                                
!c      IF(IPR.GE.3)  WRITE( 6,13 )     AICMIN                            
!x      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,11 )     I1
!x      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,13 )     AICMIN
do 330  j0=1,imp1
j1 = j0 - 1
if( j1 .eq. 0 .and. i1.ge. ii )  go to 330
mdel = (j1-1)*id + i1 + ksw
if( j1 .eq. 0 )     mdel = md + i1
do 310  i=1,md2
if( jnd(i) .eq. mdel )  go to 320
310 continue
go to 330
320 k = k0 - 1
if( i .gt. k0 )   go to 400
!c      CALL  DELETE( Y,D,IND,JND,KK,K0,MDEL,MJ4 )                        
call  delete( y,ind,jnd,kk,k0,mdel,mj4 )
!                                                                       
call  aiccom( y,n,k,kk,mj4,osd,oaic )
!                                                                       
!c      IF(IPR.GE.3)  WRITE( 6,7 )   MDEL,OAIC,MDEL,J1,I1                 
if( oaic .ge. aicmin )  go to 340
!                                                                       
!c      IF( IPR .GE. 3 )     WRITE( 6,8 )                                 
!x      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,78 ) MDEL,OAIC,MDEL,J1,I1
aicmin = oaic
call  copy( y,kk1,0,0,mj4,mj1,x )
!                                                                       
330 continue
go to 400
!c  340 K0 = K0 + 1                                                       
340 continue
k0 = k0 + 1
!c      IF( IPR .GE. 3 )     WRITE( 6,18 )                                
!x      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,718 ) MDEL,OAIC,MDEL,J1,I1
do 350  i=1,kk
j = ind(i)
!xx  350 JND(J) = I                                                        
jnd(j) = i
350 continue
400 continue
!c      IF( IPR .GE. 3 )     WRITE( 6,16 )   AICMIN , K0                  
!x      IF( (IPR.GE.3).AND.(IFG.NE.0) )   WRITE( LU,16 )   AICMIN , K0
!                                                                       
!                                                                       
!                                                                       
!c      IF(IPR.GE.3)  WRITE( 6,15 )                                       
!c      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD )                       
!x      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,15 )
!x      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD,AAIC(II),IFG,LU )
!xx      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD,AAIC(II) )
call  srcoef( x,k0,kk,n,mj1,jnd,a,sdd,aaic(ii) )
npr(ii) = k0
do 450  kk = 1,k0
jndf(kk,ii) = jnd(kk)
af(kk,ii) = a(kk)
450 continue
!c      EX(II) = SDD / DFLOAT(N)                                          
ex(ii) = sdd / dble(n)
call  coef2( a,k0,id,ii,jnd,lmax,m,ksw,msw,mj2,b,c,e )
aicsum = aicsum + aicmin
nsw = max0( nsw,msw )
!                                                                       
!                                                                       
500 continue
!                                                                       
if( isw .eq. 1 )   return
!c      IF(IPR.GE.1)  WRITE( 6,3 )                                        
!c      CALL  MCOEF( B,C,E,EX,ID,LMAX,NSW,IPR,MJ2,MJ3 )                   
call  mcoef( bi,b,c,ei,e,ex,id,lmax,nsw,ipr,mj2,mj3 )
!c      IF(IPR.GE.1)  WRITE( 6,611 )  AICSUM                              
!                                                                       
aics = aicsum
!                                                                       
return
!cx    3 FORMAT( 1H ,132(1H-) )                                            
!c    4 FORMAT( ' *****  INITIAL ESTIMATE *****' )                        
!cx    4 FORMAT( /,' *****  INITIAL ESTIMATE *****' )
!cx    5 FORMAT( ' ADD',4X ,I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,
!cx     *' ( LAG =',I2,' , I=',I2,' ) ' )
!cx    6 FORMAT( 1H+,45X,'-----',36X,'ADDED    -----' )                    
!cx   56 FORMAT( ' ADD',4X ,I3,5X,'AIC =',F15.3,10X,'-----  REGRESSOR',I4,
!cx     *' ( LAG =',I2,' , I=',I2,' )  ADDED    -----' )
!cx    7 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,
!cx     *' ( LAG =',I2,' , I=',I2,' ) ' )
!cx    8 FORMAT( 1H+,45X,'-----',36X,'DELETED  -----' )                    
!cx   78 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,10X,'-----  REGRESSOR',I4,
!cx     *' ( LAG =',I2,' , I=',I2,' )  DELETED  -----' )
!cx    9 FORMAT( 1H+,86X,'NOT ADDED' )                                     
!cx   59 FORMAT( ' ADD',4X ,I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,
!cx     *' ( LAG =',I2,' , I=',I2,' )  NOT ADDED' )
!c   11 FORMAT( 1H ,10X,'--------  CHECK',I3,'-TH  VARIABLE  --------' )  
!c   12 FORMAT( 1H ,10X,'<<<  AIC =',F15.3,'  >>>   .....TEMPORARY MINIMUM
!cx   11 FORMAT( /10X,' --------  CHECK',I3,'-TH  VARIABLE  --------' )  
!cx   12 FORMAT( /10X,' <<<  AIC =',F15.3,'  >>>   .....TEMPORARY MINIMUM',
!cx     1' AIC.....   ( NUMBER OF PARAMETERS=',I3,' )' )
!c   13 FORMAT( 1H ,15X,'AIC =',F15.3 )                                   
!c   14 FORMAT( 1H ,'*****  SECONDARY ESTIMATE  *****' )                  
!c   15 FORMAT( 1H ,'*****  FINAL ESTIMATE  *****' )                      
!c   16 FORMAT( 1H ,11X,'<<< MINIMUM AIC =',F15.3,' >>>    ( NUMBER OF PAR
!cx   13 FORMAT( /15X,' AIC =',F15.3 )                                   
!cx   14 FORMAT( /' *****  SECONDARY ESTIMATE  *****' )                  
!cx   15 FORMAT( /' *****  FINAL ESTIMATE  *****' )                      
!cx   16 FORMAT( /11X,' <<< MINIMUM AIC =',F15.3,' >>>    ( NUMBER OF ',
!cx     *'PARAMETERS =',I3,' )' )
!cx   17 FORMAT( 1H+,22X,'<<< MAICE >>>' )                                 
!cx   18 FORMAT( 1H+,86X,'NOT DELETED' )                                   
!cx  718 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,
!cx     *' ( LAG =',I2,' , I=',I2,' )  NOT DELETED' )
!c   19 FORMAT( 1H ,10X,'--------  CHECK CONSTANT VECTOR  --------' )     
!cx   19 FORMAT( /10X,' --------  CHECK CONSTANT VECTOR  --------' )     
!cx   21 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,        
!cx     * ' (CONSTANT  VECTOR)')                                           
!cx  218 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,10X,'-----  REGRESSOR',I4,
!cx     * ' (CONSTANT  VECTOR)  DELETED  -----')
!cx 2118 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,10X,'-----  REGRESSOR',I4,
!cx     * ' (CONSTANT  VECTOR)  DELETED  -----')
!c  600 FORMAT( 1H ,'-----  X  -----' )                                   
!c  602 FORMAT( 1H ,'AICMIN =',D13.5,5X,'AIC1 =',D13.5,5X,'AIC2 =',D13.5, 
!cx  600 FORMAT( /' -----  X  -----' )                                   
!cx  602 FORMAT( /' AICMIN =',D13.5,5X,'AIC1 =',D13.5,5X,'AIC2 =',D13.5, 
!cx     1 5X,'K0 =',I5,5X,'K01 =',I5,5X,'K02 =',I5 )                       
!c  603 FORMAT( 1H ,'-----  Y  -----' )                                   
!c  605 FORMAT( 1H ,4X,'I',4X,'JND(I)',10X,'A(I)' )                       
!c  606 FORMAT( 1H ,'LMAX =',I5 )                                         
!cx  603 FORMAT( /' -----  Y  -----' )                                   
!cx  605 FORMAT( /5X,'I',4X,'JND(I)',10X,'A(I)' )                       
!cx  606 FORMAT( /' LMAX =',I5 )                                         
!cx  607 FORMAT( /I5,I10,D20.5 )                                        
!c  611 FORMAT( 1H ,'AIC =',F15.3 )                                       
!c  614 FORMAT( 1H ,'--  ADDVAR VARIABLES     MA =',I5,5X,'J1 =',I5,5X,   
!cx  611 FORMAT( /' AIC =',F15.3 )                                       
!cx  614 FORMAT( /' --  ADDVAR VARIABLES     MA =',I5,5X,'J1 =',I5,5X,   
!cx     1  'I1 =',I5 )                                                     
!c  615 FORMAT( 1H ,'--  DELETE VARIABLES     MDEL =',I5,5X,'J1 =',I5,5X, 
!cx  615 FORMAT( /' --  DELETE VARIABLES     MDEL =',I5,5X,'J1 =',I5,5X, 
!cx     1  'I1 =',I5 )                                                     
!cx  616 FORMAT( 1H ,I5,F15.3 )                                            
!cx  630 FORMAT( 1H ,40I3 )                                                
!cx  635 FORMAT( 1H ,'-----  JND  -----' )                                 
!cx  637 FORMAT( 1H ,'-----  IND  -----' )                                 
!cx  638 FORMAT( 1H ,'OAIC =',F15.3,5X,'AIC1 =',F15.3 )                    
!cx  639 FORMAT( 1H ,'OAIC =',F15.3,5X,'AIC2 =',F15.3 )                    
!cx  640 FORMAT( 1H ,4X,'I',10X,'AIC ' )                                   
!cx  641 FORMAT( 1H ,'OAIC =',F15.3,5X,'AICMIN =',F15.3 )                  
!c  642 FORMAT( 1H ,30(1H-),/,1H ,'FIRST STEP OF AIC MINIMIZATION',/,1H , 
!cx  642 FORMAT( /1H ,30(1H-),/,1H ,'FIRST STEP OF AIC MINIMIZATION',/,1H ,
!cx     1 30(1H-) )                                                        
!c  643 FORMAT( 1H ,31(1H-),/,1H ,'SECOND STEP OF AIC MINIMIZATION',/,1H ,
!cx  643 FORMAT(/1H ,31(1H-),/,1H ,'SECOND STEP OF AIC MINIMIZATION',/,1H ,
!cx     1  31(1H-) )                                                       
!c  644 FORMAT( 1H ,30(1H-),/,1H ,'FINAL STEP OF AIC MINIMIZATION',/,1H , 
!cx  644 FORMAT( /1H ,30(1H-),/,1H ,'FINAL STEP OF AIC MINIMIZATION',/,1H ,
!cx     1 30(1H-) )                                                        
!c  645 FORMAT( 1H ,23X,10(1H.),'REGRESSION MODEL FOR THE REGRESSAND  II =
!cx  645 FORMAT( /24X,10(1H.),'REGRESSION MODEL FOR THE REGRESSAND  II =
!cx     1',I2,2X,10(1H.) )                                                 
end
!
!
!c      SUBROUTINE  MBYSAR( X,D,N,M,ID,KSW,IPR,MJ1,MJ2,A,B,G,H,E,AICB,EK )
!xx      SUBROUTINE  MBYSAR( X,N,M,ID,KSW,IPR,MJ1,MJ2,SD1,AIC1,DIC1,
subroutine  mbysar( x,n,m,id,ksw,mj1,mj2,sd1,aic1,dic1,&
&aicm1,sdmin1,imin1,c,d,a,b,g,h,e,aicb,ek )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PRODUCES MULTI-VARIATE AUTOREGRESSIVE MODELS BY A 
!     BAYESIAN PROCEDURE USING THE OUTPUT OF SUBROUTINE MREDCT.         
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             COPY                                                      
!             MAICE                                                     
!             BAYSWT                                                    
!             HUSHLD                                                    
!             HUSHL1                                                    
!             MARCOF                                                    
!             MBYSPC                                                    
!             MPARCO                                                    
!             MSDCOM                                                    
!             PRINT3                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!         X:      ((M+1)*ID+KSW)*((M+1)*ID+KSW) TRIANGULAR MATRIX, OUTPU
!                 SUBROUTINE MREDCT                                     
!         D:      WORKING AREA                                          
!         N:      DATA LENGTH                                           
!         M:      MAXIMUM TIME LAG OF THE MODEL                         
!         ID:     DIMENSION OF THE OBSERVATION                          
!         KSW:    =0   CONSTANT TERM IS NOT INCLUDED AS A REGRESSOR     
!                 =1   CONSTANT TERM IS INCLUDED AS THE FIRST REGRESSOR 
!         IPR:    PRINT OUT CONTROL                                     
!         MJ1:    ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM           
!         MJ2:    ABSOLUTE DIMENSION OF A,B,G,H AND E IN THE MAIN PROGRA
!                                                                       
!       OUTPUTS:                                                        
!         A:      AR-COEFFICIENT MATRICES OF FORWARD MODEL              
!         B:      AR-COEFFICIENT MATRICES OF BACKWARD MODEL             
!         G:      PARTIAL AUTOCORRELATION COEFFICIENT MATRICES OF FORWAR
!         H:      PARTIAL AUTOCORRELATION COEFFICIENT MATRICES OF BACKWA
!         E:      INNOVATION VARIANCE-COVARIANCE MATRIX                 
!         AICB:   EQUIVALENT AIC OF THE BAYESIAN (FORWARD) MODEL        
!         EK:     EQUIVALENT NUMBER OF AUTOREGRESSIVE COEFFICIENTS      
!                                                                       
!xx      IMPLICIT REAL * 8  ( A-H , O-Z )                                  
!c      DIMENSION  X(MJ1,1) , D(1) , A(MJ2,MJ2,1) , B(MJ2,MJ2,1)          
!x      DIMENSION  X(MJ1,1), D(M), A(MJ2,MJ2,1), B(MJ2,MJ2,1)          
!x      DIMENSION  G(MJ2,MJ2,1) , H(MJ2,MJ2,1)                            
!x      DIMENSION  E(MJ2,1)                                               
!xx      DIMENSION  X(MJ1,(M+1)*ID+KSW), D(M), A(MJ2,MJ2,M), B(MJ2,MJ2,M)          
!xx      DIMENSION  G(MJ2,MJ2,M) , H(MJ2,MJ2,M)                            
!xx      DIMENSION  E(MJ2,ID)
!c      DIMENSION  AIC(51) , SD(51)                                       
!xx      DIMENSION  AIC1(M+1) , SD1(M+1), DIC1(M+1), DIC2(M+1)
!xx      DIMENSION  AIC(M+1) , SD(M+1)
!c      DIMENSION  IND(100) , JND(100) , C(100)                           
!c      DIMENSION  Y(100,10)                                              
!xx      DIMENSION  IND((M+1)*ID+KSW) , JND((M+1)*ID+KSW)
!xx      DIMENSION  C(M+1)
!xx      DIMENSION  Y((M+1)*ID+KSW,ID)                                     
integer n, m, id, ksw, mj1, mj2, imin1
real(dp) x(mj1,(m+1)*id+ksw), sd1(m+1), aic1(m+1),&
&dic1(m+1), aicm1, sdmin1, c(m+1), d(m),&
&a(mj2,mj2,m), b(mj2,mj2,m), g(mj2,mj2,m),&
&h(mj2,mj2,m), e(mj2,id), aicb, ek
! local
integer i, ii, i1, imin, isub, j, j0, j1, j2, j3, jj, k, k1, m1,&
&mj4, md, md2, mdmk, ind((m+1)*id+ksw), jnd((m+1)*id+ksw)
real(dp) dic2(m+1), aic(m+1), sd(m+1), y((m+1)*id+ksw,id),&
&yy, aicm, sdmin, osd
!
!                                                                       
!          ---------------                                              
!          INITIAL SETTING                                              
!          ---------------                                              
!                                                                       
!c      MJ4 = 100                                                         
m1 = m + 1
md = m * id + ksw
md2= m1*id + ksw
mj4 = md2
!                                                                       
!xx      DO 30  I=1,M1                                                     
!xx      SD(I)  = 1.D0                                                     
!xx   30 AIC(I) = 0.D0
sd(1:m1)  = 1.d0
aic(1:m1) = 0.d0
call  copy( x,md2,0,md2,mj1,mj1,x )
!                                                                       
!                                                                       
!c      IF(IPR.GE.3)  WRITE( 6,3 )                                        
!                                                                       
!          ------------------------------------------------------------ 
!          PARTIAL AUTOREGRESSION COEFFICIENT-MATRICES OF FORWARD MODEL 
!          ------------------------------------------------------------ 
call  mparco( x,id,m,ksw,0,mj1,mj2,g,h )
!c      IF( IPR .GE. 3 )     WRITE( 6,4 )                                 
!c      IF( IPR .GE. 3 )     CALL PRINT3( G,ID,ID,M,MJ2,MJ2,1 )           
!                                                                       
!          -----------------------------------------------              
!          AIC COMPUTATION ( FORWARD AND BACKWARD MODELS )              
!          -----------------------------------------------              
j0 = md
!c      ASSIGN 70 TO ISUB                                                 
isub =  70
!                                                                       
!xx   40 DO 60  II=1,M1                                                    
40 do 61  ii=1,m1
k = (ii-1)*id + ksw
k1= k+1
!xx      DO 50  J=1,ID                                                     
do 51  j=1,id
jj = j + j0
do 50  i=k1,md2
i1 = i - k
!xx   50 Y(I1,J) = X(I,JJ)                                                 
y(i1,j) = x(i,jj)
50 continue
51 continue
!                                                                       
mdmk = md2 - k
!c      CALL  HUSHLD( Y,D,MJ4,MDMK,ID )                                   
call  hushld( y,mj4,mdmk,id )
!                                                                       
do 60  i=1,id
yy = y(i,i)**2 / n
sd(ii)  = sd(ii) * yy
!xx   60 AIC(II) = AIC(II) + N*DLOG( YY ) + 2.D0*(K+1)                     
aic(ii) = aic(ii) + n*dlog( yy ) + 2.d0*(k+1)
60 continue
61 continue
!c      IF( IPR .GE. 2 )     WRITE( 6,615 )                               
!c      GO TO ISUB, ( 70,100 )                                            
if( isub .eq. 70 ) go to 70
if( isub .eq. 100 ) go to 100
!                                                                       
!c   70 CALL  MAICE( AIC,SD,M,IPR-1,AICM,SDMIN,IMIN )                     
!xx   70 CALL  MAICE( AIC,SD,M,IPR-1,AICM,SDMIN,IMIN,DIC1 )
70 call  maice( aic,sd,m,aicm,sdmin,imin,dic1 )
do 71  i=1,m1
aic1(i) = aic(i)
sd1(i) = sd(i)
71 continue
aicm1 = aicm
sdmin1 = sdmin
imin1 = imin
!          -------------------------------------------------------------
!          PARTIAL AUTOREGRESSION COEFFICIENT-MATRICES OF BACKWARD MODEL
!          -------------------------------------------------------------
do 80  i=1,md2
!xx   80 IND(I) = I                                                        
ind(i) = i
80 continue
jnd(1) = 1
j2 = ksw
!xx      DO 85  JJ=2,M                                                     
do 86  jj=2,m
j = (m-jj)*id
do 85  i=1,id
j2 = j2 + 1
j1 = j +i
!xx   85 JND(J2) = J1                                                      
jnd(j2) = j1
85 continue
86 continue
do 90  i=1,id
j2 = j2 + 1
j1 = md + i
jnd(j2) = j1
j3 = j2 + id
j1 = j1 - id
!xx   90 JND(J3) = J1                                                      
jnd(j3) = j1
90 continue
!                                                                       
!c      CALL  HUSHL1( X,D,MJ1,MD2,MD2,1,IND,JND )                         
call  hushl1( x,mj1,md2,md2,1,ind,jnd )
call  mparco( x,id,m,ksw,1,mj1,mj2,g,h )
!c      IF( IPR .GE. 3 )      WRITE( 6,5 )                                
!c      IF( IPR .GE. 3 )      WRITE( 6,9 )                                
!c      IF( IPR .GE. 3 )      CALL PRINT3( H,ID,ID,M,MJ2,MJ2,1 )          
!                                                                       
!          ---------------                                              
!          AIC COMPUTATION                                              
!          ---------------                                              
j0 = md - id
!c      ASSIGN 100 TO ISUB                                                
isub = 100
go to 40
!                                                                       
100 do 110  i=1,m1
sd(i) = dsqrt( sd(i) )
!xx  110 AIC(I) = AIC(I) / 2                                               
aic(i) = aic(i) / 2
110 continue
!                                                                       
!c      IF( IPR .GE. 3 )     WRITE( 6,616 )                               
!c      CALL  MAICE( AIC,SD,M,IPR-2,AICM,SDMIN,IMIN )                     
!xx      CALL  MAICE( AIC,SD,M,IPR-2,AICM,SDMIN,IMIN,DIC2 )
call  maice( aic,sd,m,aicm,sdmin,imin,dic2 )
!                                                                       
!          ----------------------------                                 
!          BAYESIAN WEIGHTS COMPUTATION                                 
!          ----------------------------                                 
!                                                                       
call  bayswt( aic,aicm,m,0,c )
!                                                                       
!          -------------------------------------------------            
!          PARTIAL AUTOREGRESSION MATRICES OF BAYESIAN MODEL            
!          -------------------------------------------------            
!                                                                       
!xx      CALL  MBYSPC( G,H,C,D,M,ID,IPR,MJ2 )                              
call  mbyspc( g,h,c,d,m,id,mj2 )
!c      IF( IPR .GE. 1 )     WRITE( 6,11 )                                
!c      IF( IPR .LT. 2 )     GO TO 115                                    
!c      WRITE( 6,4 )                                                      
!c      CALL  PRINT3( G,ID,ID,M,MJ2,MJ2,1 )                               
!c      WRITE( 6,5 )                                                      
!c      CALL  PRINT3( H,ID,ID,M,MJ2,MJ2,1 )                               
!                                                                       
!          -----------------------------------------------------        
!          AUTOREGRESSION COEFFICIENT MATRICES OF BAYESIAN MODEL        
!          -----------------------------------------------------        
!                                                                       
!xx  115 CALL  MARCOF( G,H,ID,M,MJ2,A,B )                                  
call  marcof( g,h,id,m,mj2,a,b )
!c      IF( IPR .GE. 1 )   WRITE( 6,6 )                                   
!c      IF( IPR .GE. 1 )   CALL  PRINT3( A,ID,ID,M,MJ2,MJ2,1 )            
!c      IF( IPR .GE. 3 )   WRITE( 6,7 )                                   
!c      IF( IPR .GE. 3 )   CALL  PRINT3( B,ID,ID,M,MJ2,MJ2,1 )            
!          -------------------------                                    
!          AIC OF THE BAYESIAN MODEL                                    
!          -------------------------                                    
ek = 0.d0
do 120  i=1,m
!c  120 EK = EK + C(I)**2                                                 
!xx  120 EK = EK + D(I)**2                                                 
ek = ek + d(i)**2
120 continue
ek = ek*(id**2)
call  copy( x,md2,md2,0,mj1,mj1,x )
!c      CALL  MSDCOM( X,A,Y,D,N,M,ID,KSW,IPR,MJ1,MJ2,MJ4,E,OSD )          
!xx      CALL  MSDCOM( X,A,N,M,ID,KSW,IPR,MJ1,E,OSD )          
call  msdcom( x,a,n,m,id,ksw,mj1,e,osd )
!                                                                       
aicb = n*dlog( osd ) + 2.d0*ek + 2.d0*ksw*id + id*(id+1)
!c      IF( IPR .GE. 1 )     WRITE( 6,614 )   AICB                        
!                                                                       
return
!cx    3 FORMAT( 1H ,132(1H-) )                                            
!cx    4 FORMAT( /,1H ,'PARTIAL AUTOREGRESSION COEFFICIENTS  ( FORWARD MODEL
!cx     1)' )                                                              
!cx    5 FORMAT( /,1H ,'PARTIAL AUTOREGRESSION COEFFICIENTS ( BACKWARD MODEL
!cx     1)' )                                                              
!cx    6 FORMAT( /,1H ,'AR-COEFFICIENT MATRICES  ( FORWARD MODEL )' )      
!cx    7 FORMAT( /,1H ,'AR-COEFFICIENT MATRICES  ( BACKWARD MODEL )' )     
!cx    8 FORMAT( 1H ,10X,'SD  = RESIDUAL VARIANCE',15X,'=',D23.12,/,11X,'EK
!cx     1  = EQUIVALENT NUMBER OF PARAMETERS =',F15.3,/,11X,'AIC =  N*LAG(S
!cx     2D) + 2*EK',15X,'=',F15.3 )                                        
!cx    9 FORMAT( /,1H ,23(1H-),/,' LEAST SQUARES ESTIMATES',/,1H ,23(1H-) )
!cx   11 FORMAT( /,1H ,18(1H-),/,' BAYESIAN ESTIMATES',/,1H ,18(1H-) )     
!cx  610 FORMAT( 1H ,10D13.5 )                                             
!cx  612 FORMAT( 1H ,'-----  LEAST SQUARE ESTIMATES OF FULL ORDER MODEL  --
!cx     1---' )                                                            
!cx  613 FORMAT( 1H ,'-----  BAYESIAN ESTIMATES OF REGRESSION COEFFICIENTS 
!cx     * -----' )                                                         
!cx  614 FORMAT( 1H ,10X,'AIC = ',F15.3)                                   
!cx  615 FORMAT( /,1H ,'-----  FORWARD MODELS  -----' )                    
!cx  616 FORMAT( /,1H ,'-----  AVERAGE OF FORWARD AND BACKWARD MODELS  
!cx     :-----' )                                                          
!cx  642 FORMAT( 1H ,30(1H-),/,1H ,'FIRST STEP OF AIC MINIMIZATION',/,1H , 
!cx     1 30(1H-) )                                                        
!cx  645 FORMAT( 1H ,23X,10(1H.),2X,'REGRESSION MODEL FOR THE REGRESSAND  I
!cx     1I =',I2,2X,10(1H.) )                                              
end
!
!
!xx      SUBROUTINE  MBYSPC( G,H,C,D,M,ID,IPR,MJ2 )                        
subroutine  mbyspc( g,h,c,d,m,id,mj2 )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PRODUCES PARTIAL AUTOREGRESSION COEFFICIENTS G(I),
!     (I=1,K) OF THE MULTI-VARIATE AUTOREGRESSIVE MODEL.                
!                                                                       
!       INPUTS:                                                         
!          G:      LEAST SQUARES ESTIMATES OF "PARCOR'S" (FORWARD MODEL)
!          H:      LEAST SQUARES ESTIMATES OF "PARCOR'S" (BACKWARD MODEL
!          C(I+1): BAYESIAN WEIGHT OF EACH ORDER  (I=0,...,M)           
!          M:      MAXIMUM TIME LAG OF THE MODEL                        
!          ID:     DIMENSION OF THE PROCESS                             
!          IPR:    PRINT OUT CONTROL                                    
!          MJ2:    ABSOLUTE DIMENSION OF G AND H                        
!                                                                       
!       OUTPUTS:                                                        
!          G:      BAYESIAN ESTIMATES OF "PARCOR'S" (FORWARD MODEL)     
!          H:      BAYESIAN ESTIMATES OF "PARCOR'S" (BACKWARD MODEL)    
!          D:      INTEGRATED BAYESIAN WEIGHT EACH ORDER (I=1,...,M)    
!                                                                       
!xx      IMPLICIT  REAL * 8 ( A-H , O-Z )                                  
!x      DIMENSION  G(MJ2,MJ2,1) , H(MJ2,MJ2,1) , C(1) , D(1)              
!xx      DIMENSION  G(MJ2,MJ2,M) , H(MJ2,MJ2,M) , C(M+1) , D(M)
integer m, id, mj2
real(dp) g(mj2,mj2,m), h(mj2,mj2,m), c(m+1), d(m)
! local
integer i, ii, j, m1
!                                                                       
!          INTEGRATED BAYESIAN WEIGHT                                   
!                                                                       
m1 = m+1
d(m) = c(m1)
do 10 i=2,m
j = m1 - i
!xx   10 D(J) = D(J+1) + C(J+1)                                            
d(j) = d(j+1) + c(j+1)
10 continue
!c      IF( IPR .LE. 1 )   GO TO 20                                       
!c      WRITE( 6,6 )                                                      
!c      DO 15  I=1,M1                                                     
!c       IM1 = I - 1                                                      
!c      IF( I .EQ. 1 )   WRITE( 6,7 )   IM1 , C(I)                        
!c      IF( I .NE. 1 )   WRITE( 6,7 )   IM1 , C(I) , D(IM1)               
!c   15 CONTINUE                                                          
!                                                                       
!c   20 DO 30  I=1,M                                                      
!c   30 C(I) = D(I)                                                       
!                                                                       
!          PARTIAL CORRELATION                                          
!                                                                       
!xx      DO 40  II=1,M                                                     
!xx      DO 40  J=1,ID                                                     
do 42  ii=1,m
do 41  j=1,id
do 40  i=1,id
g(i,j,ii) = g(i,j,ii)*d(ii)
!xx   40 H(I,J,II) = H(I,J,II)*D(II)                                       
h(i,j,ii) = h(i,j,ii)*d(ii)
40 continue
41 continue
42 continue
!                                                                       
!                                                                       
return
!xx    6 FORMAT( 1H ,4X,'M',14X,'BAYESIAN WEIGHTS',8X,'INTEGRATES BAYESIAN 
!xx     1WEIGHTS' )                                                        
!xx    7 FORMAT( 1H ,I5,2D30.7 )                                           
end
!
!
!c      SUBROUTINE  MCOEF( B,C,E,EX,ID,LMAX,KSW,IPR,MJ2,MJ3 )             
subroutine  mcoef( bi,b,c,ei,e,ex,id,lmax,ksw,ipr,mj2,mj3 )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE COMPUTES AND PRINTS OUT THE COEFFICIENT MATRICES O
!     MULTI-VARIATE AUTOREGRESSIVE MODEL FROM THE COEFFICIENT MATRICES O
!     THE MODEL WITH INSTANTANEOUS RESPONSE                             
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             INVDET                                                    
!             PRINT3                                                    
!             TRIINV                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          B:     REGRESSION COEFFICIENTS MATRIX                        
!          C:     CONSTANT VECTOR                                       
!          E:     COEFFICIENT OF QUICK RESPONSE                         
!          EX:    RESIDUAL VARIANCES OF ORTHOGONALIZED MODEL            
!          ID:    DIMENSION OF THE OBSERVATION                          
!          LMAX:  HIGHEST ORDER OF THE AR-MODEL                         
!          KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR    
!                 =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR
!          IPR:   PRINT OUT CONTROL                                     
!          MJ2:   ABSOLUTE DIMENSION OF B AND E                         
!          MJ3:   ABSOLUTE DIMENSION OF B                               
!       OUTPUTS:                                                        
!          B:     AR-COEFFICIENT MATRIX                                 
!          C:     CONSTANT VECTOR                                       
!          E:     INNOVATION COVARIANCE MATRIX                          
!                                                                       
!xx      IMPLICIT  REAL*8 ( A-H,O-Z )                                      
!x      DIMENSION  B( MJ2,MJ2,MJ3 ) , E( MJ2,1 ) , EX( 1 ) , C( 1 )       
!c      DIMENSION  C1(24) , EE(24,24)                                     
!xx      DIMENSION  B( MJ2,MJ2,MJ3 ) , E( MJ2,ID ) , EX( ID ) , C( ID )       
!xx      DIMENSION  C1(ID) , EE(ID,ID)
integer id, lmax, ksw, ipr, mj2, mj3
real(dp) bi(id,id,lmax), b(mj2,mj2,mj3), c(id), ei(id,id),&
&e(mj2,id), ex( id )
! local
integer i, ii, j, jj, mj5
real(dp) c1(id), ee(id,id), sum, edet
!
!    INPUT  E ---> EI
!    INPUT  B ---> BI
!xx      DIMENSION  EI(ID,ID), BI(ID,ID,LMAX)
!
!c      MJ5 = 24                                                          
mj5 = id
!                                                                       
if(ipr.le.1)  go to 300
!c      WRITE( 6,603 )                                                    
!xx      DO 310  I=1,ID                                                    
do 311  i=1,id
do 310  j=1,id
!xx  310 EI(I,J) = E(I,J)
ei(i,j) = e(i,j)
310 continue
311 continue
!c  310 WRITE( 6,610 )  (E(I,J),J=1,ID)                                   
!c      WRITE( 6,604 )                                                    
do  330     ii=1,lmax
!xx      DO  320     I=1,ID
do  321     i=1,id
do  320     j=1,id
!xx  320 BI(I,J,II) = B(I,J,II)
bi(i,j,ii) = b(i,j,ii)
320 continue
321 continue
!c  320 WRITE( 6,610 )     (B(I,J,II),J=1,ID)                             
!c  330 WRITE( 6,609 )                                                    
330 continue
!xx      IF( KSW .EQ. 0 )     GO TO 300                                    
!c      WRITE( 6,601 )                                                    
!c      WRITE( 6,610 )     (C(I),I=1,ID)                                  
300 continue
!                                                                       
if( ksw .ne. 1 )     go to 375
!xx      DO  335     I=1,ID                                                
do  336     i=1,id
do  335     j=1,id
sum = e(i,j)
do  325     ii=1,lmax
!xx  325 SUM = SUM - B(I,J,II)                                             
sum = sum - b(i,j,ii)
325 continue
!xx  335 EE(I,J) = SUM
ee(i,j) = sum
335 continue
336 continue
call  invdet( ee,edet,id,mj5 )
do  355     i=1,id
sum = 0.d0
do  345     j=1,id
!xx  345 SUM = SUM + EE(I,J)*C(J)                                          
sum = sum + ee(i,j)*c(j)
345 continue
!xx  355 C1(I) = SUM                                                       
c1(i) = sum
355 continue
do  365     i=1,id
!xx  365 C(I) = C1(I)                                                      
c(i) = c1(i)
365 continue
!                                                                       
375 continue
call  triinv( e,id,mj2,mj5,ee )
!                                                                       
!xx      DO  360     II=1,LMAX                                             
!xx      DO  350     I=1,ID
do  362     ii=1,lmax
do  351     i=1,id
do  350     j=1,id
sum = 0.d0
do  340     jj=1,i
!xx  340 SUM = SUM + EE(I,JJ)*B(JJ,J,II)
sum = sum + ee(i,jj)*b(jj,j,ii)
340 continue
!xx  350 E(I,J) = SUM
e(i,j) = sum
350 continue
351 continue
!xx      DO  360     I=1,ID                                                
do  361     i=1,id
do  360     j=1,id
!xx  360 B(I,J,II) = E(I,J)                                                
b(i,j,ii) = e(i,j)
360 continue
361 continue
362 continue
do  370     i=1,id
sum = 0.d0
do  371     j=1,id
!xx  371 SUM = SUM + EE(I,J)*C(J)                                          
sum = sum + ee(i,j)*c(j)
371 continue
!xx  370 C1(I) = SUM                                                       
c1(i) = sum
370 continue
do  372     i=1,id
!xx  372 C(I) = C1(I)
c(i) = c1(i)
372 continue
!xx      DO  380     I=1,ID                                                
do  381     i=1,id
do  380     j=1,i
sum = 0.d0
do  385     ii=1,j
!xx  385 SUM = SUM + EE(I,II) * EE(J,II)*EX(II)                            
sum = sum + ee(i,ii) * ee(j,ii)*ex(ii)
385 continue
e(i,j) = sum
!xx  380 E(J,I) = SUM                                                      
e(j,i) = sum
380 continue
381 continue
!c      IF(IPR.EQ.0)   RETURN                                             
!c      WRITE( 6,612 )                                                    
!c      DO  390     I=1,ID                                                
!c  390 WRITE( 6,611 )     (E(I,J),J=1,ID)                                
!c      WRITE( 6,613 )                                                    
!c      CALL  PRINT3( B,ID,ID,LMAX,MJ2,MJ2,0 )                            
!c      IF(KSW.NE.1)  RETURN                                              
!c      WRITE( 6,602 )                                                    
!c      WRITE( 6,610 )     (C(I),I=1,ID)                                  
!                                                                       
return
!cx  601 FORMAT( 1H ,'-----  C  -----' )                                   
!cx  602 FORMAT( 1H ,'-----  CONSTANT VECTOR  -----' )                     
!cx  603 FORMAT( 1H ,'-----  INSTANTANEOUS RESPONSE  -----' )              
!cx  604 FORMAT( 1H ,'-----  B  -----' )                                   
!cx  609 FORMAT( 1H  )                                                     
!cx  610 FORMAT( 1H ,8F15.8 )                                              
!cx  611 FORMAT( 1H ,8D15.7 )                                              
!cx  612 FORMAT( 1H ,'*  INNOVATION VARIANCE MATRIX  *' )                  
!cx  613 FORMAT( 1H ,'*  AR-COEFFICIENTS  *' )                             
end
!
!
subroutine  mparco( x,id,m,ksw,ifg,mj1,mj3,g,h )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PRODUCES LEAST SQUARES ESTIMATES OF PARTIAL       
!     AUTOREGRESSION COEFFICIENT MATRICES OF FORWARD AND BACKWARD       
!     MULTI-DIMENSIONAL AUTOREGRESSIVE MODEL                            
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
!             SOLVE                                                     
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          X:      ((M+1)*ID+KSW)*((M+1)*ID+KSW) MATRIX                 
!          ID:     DIMENSION OF OBSERVATION                             
!          M:      MAXIMUM TIME LAG OF THE MODEL                        
!          KSW:    =0  CONSTANT TERM IS NOT INCLUDED AS A REGRESSOR     
!                  =1  CONSTANT TERM IS INCLUDED AS THE FIRST REGRESSOR 
!          IFG:    =0   TO REQUEST THE COEFFICIENTS OF FORWARD MODEL    
!                  =1   TO REQUEST THE COEFFICIENTS OF BACKWARD MODEL   
!          MJ1:    ABSOLUTE DIMENSION OF X                              
!          MJ3:    ABSOLUTE DIMENSION OF G AND H                        
!                                                                       
!       OUTPUTS:                                                        
!          G:      AR-COEFFICIENT MATRICES OF FORWARD MODEL             
!          H:      AR-COEFFICIENT MATRICES OF BACKWARD MODEL            
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  X(MJ1,1) , G(MJ3,MJ3,1) , H(MJ3,MJ3,1)                 
!c      DIMENSION  C(10,10) , R(10,10)                                    
!cx      DIMENSION  X(MJ1,(M+1)*ID+KSW) , G(MJ3,MJ3,M) , H(MJ3,MJ3,M)                 
!cx      DIMENSION  C(ID,ID) , R(ID,ID)
integer id, m, ksw, ifg, mj1, mj3
real(dp) x(mj1,(m+1)*id+ksw), g(mj3,mj3,m), h(mj3,mj3,m)
! local
integer i, ii, i0, i1, j, js0, j0, j1, j2, md, mj2, mm1
real(dp) c(id,id), r(id,id)
!                                                                       
!c      MJ2 = 10                                                          
mj2 = id
if( ifg .ne. 0 )     go to 30
md = id*m + ksw
do  20     ii=1,m
i0 = (ii-1)*id + ksw
!xx      DO  10     J=1,ID                                                 
do  11     j=1,id
j1 = j + i0
j2 = j + md
do  10     i=1,id
i1 = i + i0
c(i,j) = x(i1,j1)
!xx   10 R(I,J) = X(I1,J2)
r(i,j) = x(i1,j2)
10 continue
11 continue
!xx   20 CALL  SOLVE( C,R,ID,II,MJ2,MJ3,G )                                
call  solve( c,r,id,ii,mj2,mj3,g )
20 continue
go to 60
!                                                                       
30 mm1 = m - 1
js0 = mm1*id + ksw
do  50     ii=1,m
i0 = (ii-1)*id + ksw
j0 = (mm1-ii)*id + ksw
if( ii .eq. m )     j0 = m*id + ksw
!xx      DO  40     J=1,ID
do  41     j=1,id
j1 = j0 + j
j2 = js0 + j
do  40     i=1,id
i1 = i0 + i
c(i,j) = x(i1,j1)
!xx   40 R(I,J) = X(I1,J2)                                                 
r(i,j) = x(i1,j2)
40 continue
41 continue
!xx   50 CALL  SOLVE( C,R,ID,II,MJ2,MJ3,H )                                
call  solve( c,r,id,ii,mj2,mj3,h )
50 continue
!                                                                       
60 return
end
!
!
!c      SUBROUTINE  MRDATA( MT,MJ,Z,N,ID )                                
subroutine mrdata( zs,z,n,id,c,zmean,zvari )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!         +-----------------------------------------+                   
!         ! ORIGINAL DATA LOADING AND MEAN DELETION !                   
!         +-----------------------------------------+                   
!                                                                       
!     THIS SUBROUTINE IS USED FOR THE LOADING OF ORIGINAL DATA AND      
!     DELETION OF THE MEAN VALUES.  THE DATA IS LOADED THROUGH THE DEVIC
!     SPECIFIED BY MT.  EACH DATA SET IS COMPOSED OF TITLE, DATA LENGTH,
!     CALIBRATIONS OF DATA, DATA FORMAT AND THE ORIGINAL DATA.          
!                                                                       
!       INPUTS:                                                         
!          MT:      INPUT DEVICE SPECIFICATION                          
!          MJ:      ABSOLUTE DIMENSION OF Z                             
!                                                                       
!          TITLE:   SPECIFICATION OF DATA                               
!          N:       DATA LENGTH                                         
!          ID:      DIMENSION OF OBSERVATION                            
!          IFM:     CONTROL FOR INPUT                                   
!          FORM:    INPUT DATA FORMAT SPECIFICATION                     
!          C(I):    CALIBRATION OF CHANNEL I (I=1,ID)                   
!          Z:       ORIGINAL DATA; Z(K,I) (K=1,N) REPRESENTS THE I-TH CH
!                   RECORD                                              
!                                                                       
!       OUTPUTS:                                                        
!          Z:       ORIGINAL DATA ( MEAN DELETED )                      
!          N:       DATA LENGTH                                         
!          ID:      DIMENSION OF OBSERVATION                            
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4  Z(MJ,1)                                                 
!c      DIMENSION  Z(MJ,1)
!c      DIMENSION  C(24)                                                  
!xx      DIMENSION  ZS(N,ID), Z(N,ID), C(ID)
!xx      DIMENSION  ZMEAN(ID), ZVARI(ID)
integer n, id
real(dp) zs(n,id), z(n,id), c(id), zmean(id), zvari(id)
! local
integer i, j
real(dp) cc, sum
!c      REAL * 4  FORM(20) , TITLE(20)                                    
!                                                                       
!                                                                       
!c      READ( MT,2 )     TITLE                                            
!c      READ( MT,1 )     N , ID , IFM                                     
!c      READ( MT,2 )     FORM                                             
!c      READ( MT,3 )     (C(I),I=1,ID)                                    
!                                                                       
!c      GO TO ( 10,20,30,40,50,60 ) , IFM                                 
!                                                                       
!c   10 DO  15    I=1,N                                                   
!c   15 READ( MT,FORM )     (Z(I,J),J=1,ID)                               
!c      GO TO 70                                                          
!c   20 DO  25     J=1,ID                                                 
!c   25 READ( MT,FORM )     (Z(I,J),I=1,N)                                
!c      GO TO 70                                                          
!c   30 READ( MT,FORM )     ((Z(I,J),J=1,ID),I=1,N)                       
!c      GO TO 70                                                          
!c   40 READ( MT,FORM )     ((Z(I,J),I=1,N),J=1,ID)                       
!c      GO TO 70                                                          
!c   50 READ( MT )     ((Z(I,J),J=1,ID),I=1,N)                            
!c      GO TO 70                                                          
!c   60 READ( MT )     ((Z(I,J),I=1,N),J=1,ID)                            
!c   70 CONTINUE                                                          
!                                                                       
!                                                                       
!       ORIGINAL DATA PRINT OUT                                         
!                                                                       
!c      WRITE( 6,601 )     TITLE                                          
!c      WRITE( 6,602 )     N , ID , FORM                                  
!c      WRITE( 6,603 )                                                    
!xx      DO  75     J=1,ID                                                 
do  76     j=1,id
!c      WRITE( 6,604 )     J                                              
!c   75 WRITE( 6,610 )     (Z(I,J),I=1,N)                                 
do  75     i=1,n
!xx   75 Z(I,J) = ZS(I,J)
z(i,j) = zs(i,j)
75 continue
76 continue
!                                                                       
!xx      DO  85     J=1,ID                                                 
do  86     j=1,id
cc = c(j)
do  85     i=1,n
!xx   85 Z(I,J) = Z(I,J) * CC                                              
z(i,j) = z(i,j) * cc
85 continue
86 continue
!                                                                       
!                                                                       
!c      WRITE( 6,606 )                                                    
do  120     j=1,id
!                                                                       
!     MEAN DELETION   AND  VARIANCE COMPUTATION                         
!                                                                       
sum = 0.d0
do  90     i=1,n
!xx   90 SUM = SUM + Z(I,J)                                                
sum = sum + z(i,j)
90 continue
!c      ZMEAN = SUM / DFLOAT(N)                                           
zmean(j) = sum / dble(n)
do  100     i=1,n
!c  100 Z(I,J) = Z(I,J) - ZMEAN                                           
!xx  100 Z(I,J) = Z(I,J) - ZMEAN(J)
z(i,j) = z(i,j) - zmean(j)
100 continue
!                                                                       
sum = 0.d0
do  110     i=1,n
!xx  110 SUM = SUM + Z(I,J)*Z(I,J)                                         
sum = sum + z(i,j)*z(i,j)
110 continue
!c      ZVARI = SUM / DFLOAT(N)                                           
zvari(j) = sum / dble(n)
!                                                                       
!         MEAN AND VARIANCE PRINT OUT                                   
!                                                                       
!c  120 WRITE( 6,605 )     J , ZMEAN , ZVARI                              
120 continue
!                                                                       
!                                                                       
return
!                                                                       
!cx    1 FORMAT( 16I5 )                                                    
!cx    2 FORMAT( 20A4 )                                                    
!cx    3 FORMAT( 8F10.0 )                                                  
!cx  601 FORMAT( 1H ,'TITLE  ---  ',20A4 )                                 
!cx  602 FORMAT( 1H ,'N =',I5,5X,'ID =',I5,5X,'FORMAT =',20A4 )            
!cx  603 FORMAT( 1H ,'** ORIGINAL DATA **' )                               
!cx  604 FORMAT( 1H ,I5,'-CHANNEL' )                                       
!cx  605 FORMAT( 1H ,I6,2D15.8 )                                           
!cx  606 FORMAT( 1H ,5X,'I',10X,'MEAN',7X,'VARIANCE' )                     
!cx  610 FORMAT( 1H ,10D13.5 )                                             
end
!
!
!c      SUBROUTINE  MREDCT( Z,D,NMK,N0,LAG,ID,MJ,MJ1,KSW,X )              
subroutine  mredct( z,nmk,n0,lag,id,mj,mj1,ksw,x )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!         +-------------------------+                                   
!         ! HOUSEHOLDER'S REDUCTION !                                   
!         +-------------------------+                                   
!                                                                       
!     THIS SUBROUTINE FIRST SETS UP DATA MATRIX X BY AUGMENTING         
!     SUCCESSIVELY SHIFTED ORIGINAL DATA MATRIX Z AND THEN TRANSFORMS X 
!     INTO TRIANGULAR FORM BY HOUSEHOLDER TRANSFORMATION.               
!                                                                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             HUSHLD                                                    
!             MSETX1                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          Z:      ORIGINAL DATA MATRIX                                 
!          D:      WORKING AREA                                         
!          NMK:    DIMENSION OF THE VECTOR OF REGRESSAND (Z(N0+K-KSW+1),
!                  ...,Z(N0+K-KSW+NMK) (J=1,ID)                         
!          N0:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
!          LAG:    HIGHEST TIME LAG OF THE MODEL                        
!          ID:     DIMENSION OF OBSERVATIONS                            
!          MJ:     ABSOLUTE DIMENSION OF Z                              
!          MJ1:    ABSOLUTE DIMENSION OF X                              
!          KSW:    =0   CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR  
!                  =1   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESS
!                                                                       
!         OUTPUT:                                                       
!          X:      REDUCED MATRIX ( UPPER TRIANGULAR FORM )             
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  X(MJ1,1) , D(1)                                        
!C      REAL * 4  Z(MJ,1)                                                 
!x      DIMENSION  Z(MJ,1), X(MJ1,1)
!xx      DIMENSION  Z(MJ,ID), X(MJ1,(LAG+1)*ID+KSW)
integer nmk, n0, lag, id, mj, mj1, ksw
real(dp) z(mj,id), x(mj1,(lag+1)*id+ksw)
! local
integer k1, kd1, l, lk, n1, n2
!                                                                       
l = min0( nmk,mj1 )
k1 = lag + 1
kd1 = k1*id + ksw
n1 = l
!                                                                       
!          +-----------------+                                        +-
!          ! MATRIX X SET UP !                                        ! 
!          +-----------------+                                        +-
!                                                                       
call  msetx1( z,n0,l,lag,id,mj,mj1,0,ksw,x )
!                                                                       
!          +----------------------------+                             +-
!          ! HOUSEHOLDER TRANSFORMATION !                             ! 
!          +----------------------------+                             +-
!                                                                       
!c      CALL  HUSHLD( X,D,MJ1,L,KD1 )                                     
call  hushld( x,mj1,l,kd1 )
!                                                                       
if( n1 .ge. nmk )     return
!                                                           +------->>--
!                                                           !           
100 l = min0( nmk-n1,mj1-kd1 )
!                                                           !           
lk = l + kd1
n2 = n0 + n1
!                                                           !           
!          +-----------------+                              !         +-
!          ! MATRIX X SET UP !                              !         ! 
!          +-----------------+                              !         +-
!                                                           !           
call  msetx1( z,n2,l,lag,id,mj,mj1,1,ksw,x )
!                                                           !           
!          +----------------------------+                   !         +-
!          ! HOUSEHOLDER TRANSFORMATION !                   !         ! 
!          +----------------------------+                   !         +-
!                                                           !           
!c      CALL  HUSHLD( X,D,MJ1,LK,KD1 )                                    
call  hushld( x,mj1,lk,kd1 )
!                                                           !         +-
!                                                           +<-- NO --!N
!                                                                     +-
n1 = n1 + l
if( n1 .lt. nmk )     go to 100
!                                                                       
!                                                                     +-
!                                                                     ! 
!                                                                     +-
return
!                                                                       
end
!
!
!c      SUBROUTINE  MSDCOM( X,A,Y,D,N,M,ID,KSW,IPR,MJ,MJ2,MJ4,E,SD )      
!xx      SUBROUTINE  MSDCOM( X,A,N,M,ID,KSW,IPR,MJ,E,SD )      
subroutine  msdcom( x,a,n,m,id,ksw,mj,e,sd )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PRODUCES THE ONE-STEP AHEAD PREDICTION ERROR VARIA
!     MATRIX AND ITS DETERMINANT FOR THE MULTI-VARIATE AUTOREGRESSIVE MO
!     DEFINED BY THE AR-COEFFICIENT MATRICES GIVEN BY A.                
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
!             HUSHLD                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          X:    ((M+1)*ID+KSW)*((M+1)*ID+KSW) TRIANGULAR MATRIX        
!          A:    AR-COEFFICIENT MATRICES                                
!          Y:    WORKING AREA                                           
!          D:    WORKING AREA                                           
!          N:    DATA LENGTH                                            
!          M:    ORDER OF THE AR-MODEL                                  
!          ID:   DIMENSION OF THE PROCESS                               
!          KSW:  =0   CONSTANT TERM IS NOT INCLUDED AS A REGRESSOR      
!                =1   CONSTANT TERM IS INCLUDED AS THE FIRST REGRESSOR  
!          IPR:  PRINT OUT CONTROL                                      
!          MJ:   ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM            
!          MJ2:  ABSOLUTE DIMENSION OF A AND E IN THE MAIN PROGRAM      
!          MJ4:  ABSOLUTE DIMENSION OF Y IN THE MAIN PROGRAM            
!                                                                       
!       OUTPUTS:                                                        
!          E:    ONE STEP AHEAD PREDICTION ERROR VARIANCE MATRIX        
!          SD:   DETERMINANT OF SD                                      
!                                                                       
!xx      IMPLICIT  REAL * 8 ( A-H,O-Z )                                    
!c      DIMENSION  X(MJ,1) , A(MJ2,MJ2,1) , E(MJ2,1) , Y(MJ4,1)           
!c      DIMENSION  D(1)                                                   
!x      DIMENSION  X(MJ,1) , A(ID,ID,M) , E(ID,ID) , Y((M+1)*ID,ID)
!xx      DIMENSION X(MJ,(M+1)*ID+KSW), A(ID,ID,M), E(ID,ID), Y((M+1)*ID,ID)
integer n, m, id, ksw, mj
real(dp) x(mj,(m+1)*id+ksw), a(id,id,m), e(id,id), sd
! local
integer i, i0, i1, ii, j, j0, j1, jj, md, m0, m1d
real(dp) y((m+1)*id,id), sum
!                                                                       
md = m * id
m1d= md + id
do 30  jj=1,id
do 20  ii=1,md
i0 = ii + ksw
sum = 0.d0
do 10  j=ii,md
j0 = j + ksw
m0 = (j-1)/id + 1
j1 = j - (m0-1)*id
!xx   10 SUM = SUM + X(I0,J0)*A(JJ,J1,M0)                                  
sum = sum + x(i0,j0)*a(jj,j1,m0)
10 continue
j0 = md + ksw + jj
!xx   20 Y(II,JJ) = X(I0,J0) - SUM                                         
y(ii,jj) = x(i0,j0) - sum
20 continue
30 continue
!xx      DO 40  J=1,ID                                                     
do 41  j=1,id
j0 = md + ksw + j
do 40  i=1,id
i0 = md + ksw + i
i1 = md + i
!xx   40 Y(I1,J) = X(I0,J0)                                                
y(i1,j) = x(i0,j0)
40 continue
41 continue
!                                                                       
!c      CALL  HUSHLD( Y,D,MJ4,M1D,ID )                                    
call  hushld( y,m1d,m1d,id )
sd = 1.d0
do 50  i=1,id
!xx   50 SD = SD * (Y(I,I)**2)/N                                           
sd = sd * (y(i,i)**2)/n
50 continue
!xx      DO 70  I=1,ID                                                     
do 71  i=1,id
do 70  j=1,id
sum = 0.d0
do 60  ii=1,id
!xx   60 SUM = SUM + Y(II,I)*Y(II,J)                                       
sum = sum + y(ii,i)*y(ii,j)
60 continue
!xx   70 E(I,J) = SUM / N                                                  
e(i,j) = sum / n
70 continue
71 continue
!c      IF( IPR .LT. 1 )     RETURN                                       
!c      WRITE( 6,5 )                                                      
!c      DO 180  I=1,ID                                                    
!c  180 WRITE( 6,6 )   (E(I,J),J=1,ID)                                    
return
!cx    5 FORMAT( /,1H ,'*  INNOVATION VARIANCE MATRIX  *' )
!cx    6 FORMAT( 1H ,5D20.10 )                                             
end
!
!
subroutine  msetx1( z,n0,l,lag,id,mj,mj1,jsw,ksw,x )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          +-----------------+                                          
!          ! MATRIX X SET UP !                                          
!          +-----------------+                                          
!                                                                       
!     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA MATRIX Z(I,J) (I=
!     ,N0+K-KSW+L;J=1,...,ID) FOR AUTOREGRESSIVE MODEL FITTING.  X IS TH
!     AS INPUT TO SUBROUTINE HUSHLD.                                    
!       INPUTS:                                                         
!          Z:      ORIGINAL DATA MATRIX                                 
!          N0:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
!          L:      DIMENSION OF THE VECTOR OF NEW OBSERVATIONS          
!          LAG:    MAXIMUM TIME LAG OF THE MODEL                        
!          ID:     DIMENSION OF OBSERVATION                             
!          MJ:     ABSOLUTE DIMENSION OF Z                              
!          MJ1:    ABSOLUTE DIMENSION OF X                              
!          JSW:   =0   TO CONSTRUCT INITIAL L*(LAG+1) DATA MATRIX       
!                 =1   TO AUGMENT ORIGINAL (LAG+1)*(LAG+1) MATRIX X BY A
!                      L*(LAG+1) DATA MATRIX OF ADDTIONAL OBSERVATIONS  
!          KSW:   =0   CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR   
!                 =1   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSO
!                                                                       
!       OUTPUT:                                                         
!          X:      L*(LAG+1) MATRIX            IF  JSW = 0              
!                  (LAG+1+L)*(LAG+1) MATRIX    IF  JSW = 1              
!                                                                       
!x      REAL * 8  X(MJ1,1) , Z(MJ,1)
!xx      REAL * 8  X(MJ1,(LAG+1)*ID+KSW) , Z(MJ,ID)
integer n0, l, lag, id, mj, mj1, jsw, ksw
real(dp) z(mj,id), x(mj1,(lag+1)*id+ksw)
! local
integer i, i0, i1, i2, ii, j, j1, j2, jj, kd, kd1
!C      DIMENSION  Z(MJ,1)                                                
!                                                                       
kd = lag*id + ksw
kd1 = (lag+1)*id + ksw
i0 = 0
if( jsw .eq. 1 )     i0 = kd1
!                                                                       
!xx      DO  30     II=1,L                                                 
do  31     ii=1,l
i1 = n0 + lag + ii
i2 = i0 + ii
do  10     j=1,id
j2 = kd + j
!xx   10 X(I2,J2) = Z(I1,J)                                                
x(i2,j2) = z(i1,j)
10 continue
do  30     jj=1,lag
i1 = i1 - 1
j1 = (jj-1)*id + ksw
do  20     j=1,id
j2 = j1 + j
!xx   20 X(I2,J2) = Z(I1,J)                                                
x(i2,j2) = z(i1,j)
20 continue
30 continue
31 continue
!                                                                       
if( ksw .ne. 1 )     return
do  40     ii=1,l
i = ii + i0
!xx   40 X(I,1) = 1.D0                                                     
x(i,1) = 1.d0
40 continue
!                                                                       
return
end
!
!
!c      SUBROUTINE  NRASPE( SGME2,A,B,L,K,H,TITLE )                       
subroutine  nraspe( sgme2,a,b,l,k,h,sxx )
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES POWER SPECTRUM OF AN AR-MA PROCESS       
!     X(N)=A(1)X(N-1)+...+A(L)X(N-L)+E(N)+B(1)E(N-1)+...+B(K)E(N-K),    
!     WHERE E(N) IS A WHITE NOISE WITH ZERO MEAN AND VARIANCE EQUAL TO  
!     SGME2.  OUTPUTS PXX(I) ARE GIVEN AT FREQUENCIES I/(2*H)           
!     I=0,1,...,H.                                                      
!     REQUIRED INPUTS ARE;                                              
!     L,K,H,SGME2,(A(I),I=1,L), AND (B(I),I=1,K).                       
!        SGME2: NOISE VARIANCE                                          
!        A: AR-COEFFICIENTS                                             
!        B: MA-COEFFICIENTS                                             
!        L: ORDER OF AR                                                 
!        K: ORDER OF MA                                                 
!        N: LENGTH OF DATA                                              
!        H: NUMBER OF SEGMENTS OF FREQUENCY AXIS                        
!     0 IS ALLOWABLE AS L AND/OR K.                                     
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             FOUGER                                                    
!             SPEGRH                                                    
!       ----------------------------------------------------------------
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!C      REAL*4  SXX(501) , TITLE(1)                                       
!c      REAL*4  TITLE(1)
!xx      INTEGER H,H1                                                      
!c      DIMENSION A(501),B(501)                                           
!c      DIMENSION G(501),GR1(501),GI1(501),GR2(501),GI2(501)              
!c      DIMENSION SXX(501) , PXX(510)                                                
!xx      DIMENSION A(L),B(K)
!xx      DIMENSION G(L+K+1),GR1(H+1),GI1(H+1),GR2(H+1),GI2(H+1)
!xx      DIMENSION SXX(H+1), PXX(H+1)
integer l, k, h
real(dp) sgme2, a(l), b(k), sxx(h+1)
! local
integer i, i1, k1, l1, h1
real(dp) g(l+k+1), gr1(h+1), gi1(h+1), gr2(h+1), gi2(h+1),&
&pxx(h+1)
!xx  310 H1=H+1                                                            
h1=h+1
l1=l+1
k1=k+1
g(1)=1.0
if(l.le.0) go to 400
do 10 i=1,l
i1=i+1
!xx   10 G(I1)=-A(I)                                                       
g(i1)=-a(i)
10 continue
400 call fouger(g,l1,gr1,gi1,h1)
g(1)=1.0
if(k.le.0) go to 410
do 20 i=1,k
i1=i+1
!xx   20 G(I1)=B(I)
g(i1)=b(i)
20 continue
410 call fouger(g,k1,gr2,gi2,h1)
do 30 i=1,h1
!xx   30 PXX(I)=(GR2(I)**2+GI2(I)**2)/(GR1(I)**2+GI1(I)**2)*SGME2          
pxx(i)=(gr2(i)**2+gi2(i)**2)/(gr1(i)**2+gi1(i)**2)*sgme2
30 continue
!xx  510 DO 520 I=1,H1                                                     
do 520 i=1,h1
!xx  520 SXX(I)=DLOG10(PXX(I))                                             
sxx(i)=dlog10(pxx(i))
520 continue
!xx 1000 CONTINUE                                                          
!c      CALL  SPEGRH( SXX,TITLE )                                         
!                                                                       
return
end
!
!
subroutine  parcor( ar,k,pac )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  TRANSFORMATION FROM AR COEFFICIENTS TO PARCOR                   
!                                                                       
!       INPUTS:                                                         
!          AR:   VECTOR OF AR-COEFFICIENTS                              
!          K:    ORDER OF THE MODEL                                     
!       OUTPUT:                                                         
!          PAC:  VECTOR OF PARTIAL AUTOCORRELATIONS                     
!                                                                       
!xx      IMPLICIT  REAL * 8( A-H,O-Z )                                     
!c      DIMENSION  AR(K) , PAC(K) , W(20)                                 
!xx      DIMENSION  AR(K) , PAC(K) , W(K)
integer k
real(dp) ar(k), pac(k)
! local
integer i, i2, ii, ii1, j, jj
real(dp) w(k), s
!
do 10  i=1,k
!xx   10 PAC(I) = AR(I)                                                    
pac(i) = ar(i)
10 continue
if( k .eq. 1 )   return
do 20  jj=1,k-1
ii = k - jj
s = 1.d0 - pac(ii+1)**2
do 30  i=1,ii
j = ii - i + 1
!xx   30 W(I) = (PAC(I) + PAC(II+1)*PAC(J))/S                              
w(i) = (pac(i) + pac(ii+1)*pac(j))/s
30 continue
ii1 = ii + 1
i2 = ii1 / 2
if( mod( ii,2 ) .eq. 1 )  w(i2) = pac(i2)/(1.d0 - pac(ii1))
do 40  i=1,ii
!xx   40 PAC(I) = W(I)                                                     
pac(i) = w(i)
40 continue
20 continue
!                                                                       
return
end
!
subroutine  recoef( x,m,k,mj,a )
  use timsac_kinds, only: dp
  implicit none
!          +-------------------------------------+                      
!          ! REGRESSION COEFFICIENTS COMPUTATION !                      
!          +-------------------------------------+                      
!                                                                       
!     THIS SUBROUTINE PRODUCES REGRESSION COEFFICIENTS OF REGRESSION MOD
!     WITH M REGRESSORS FROM THE TRIANGULAR MATRIX X PREPARED BY SUBROUT
!     REDUCT.                                                           
!                                                                       
!       INPUTS:                                                         
!          X:     (K+1)*(K+1) TRIANGULAR MATRIX, WITH J-TH COLUMN REPRES
!                 J-TH REGRESSOR (J=1,K) AND (K+1)-TH REGRESSAND        
!          K:     NUMBER OF REGRESSORS                                  
!          M:     REGRESSION ON THE FIRST M REGRESSORS IS REQUESTED     
!          MJ:    ABSOLUTE DIMENSION OF X                               
!                                                                       
!       OUTPUT:                                                         
!          A(I) (I=1,M):   VECTOR OF REGRESSION COEFFICIENTS            
!                                                                       
!xx      IMPLICIT REAL * 8 (A-H,O-Z )                                      
!x      DIMENSION X(MJ,1) , A(1)                                          
!xx      DIMENSION X(MJ,K+1) , A(M)
integer m, k, mj
real(dp) x(mj,k+1), a(m)
! local
integer i, ii, i1, j, k1, mm1
real(dp) sum
!                                                                       
k1 = k + 1
a(m) = x(m, k1) / x(m,m)
if( m .eq. 1 )     return
mm1 = m - 1
do  10   ii = 1,mm1
i = m - ii
sum = x(i,k1)
i1 = i + 1
do  20   j = i1,m
!xx   20 SUM = SUM - A(J) * X(I,J)                                         
sum = sum - a(j) * x(i,j)
20 continue
!xx   10 A(I) = SUM / X(I,I)                                               
a(i) = sum / x(i,i)
10 continue
!                                                                       
return
!                                                                       
end
!
!
!c      SUBROUTINE  REDATA( X,N,MT,TITLE )                                
subroutine  redata( xs,x,n,xmean,sum )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          +---------------------------------------+                    
!          ! ORIGINAL DATA INPUT AND MEAN DELETION !                    
!          +---------------------------------------+                    
!                                                                       
!     THIS SUBROUTINE IS USED FOR THE LOADING OF ORIGINAL DATA AND DELET
!     THE MEAN VALUE.  THE DATA IS LOADED THROUGH THE DEVICE SPECIFIED B
!     EACH DATA SET IS COMPOSED OF TITLE, DATA LENGTH, DATA FORMAT AND  
!     ORIGINAL DATA.                                                    
!                                                                       
!       INPUTS:                                                         
!         MT:     INPUT DEVICE SPECIFICATION                            
!         TITLE:  TITLE OF DATA                                         
!         N:      DATA LENGTH                                           
!         DFORM:  INPUT DATA FORMAT SPECIFICATION                       
!         X(I) (I=1,N):  ORIGINAL DATA                                  
!                                                                       
!       OUTPUTS:                                                        
!         X:   ORIGINAL DATA ( MEAN DELETED )                           
!         N:   DATA LENGTH                                              
!         TITLE:  TITLE OF DATA                                         
!                                                                       
!xx      IMPLICIT REAL*8( A-H,O-Z )                                        
!C      REAL * 4     DFORM(20) , TITLE(20) , X(1)                         
!xx      DIMENSION     XS(N), X(N)
integer n
real(dp) xs(n), x(n), xmean, sum
integer i
!c	REAL * 4     DFORM(20) , TITLE(20)
!                                                                       
!       LOADING OF TITLE, DATA LENGTH, FORMAT SPECIFICATION AND DATA    
!                                                                       
!c      READ( MT,5 )     TITLE                                            
!c      READ( MT,1 )      N                                               
!c      READ( MT,5 )     DFORM                                            
!c      READ( MT,DFORM )     (X(I),I=1,N)                                 
!                                                                       
!       ORIGINAL DATA PRINT OUT                                         
!                                                                       
!c      WRITE( 6,6 )                                                      
!c      WRITE( 6,8 )     TITLE                                            
!c      WRITE( 6,9 )     N , (DFORM(I),I=1,20)                            
!c      WRITE( 6,4 )     (X(I),I=1,N)                                     
do 100 i=1,n
x(i)=xs(i)
100 continue
!                                                                  
!          MEAN DELETION                                                
!                                                                       
sum = 0.0d00
do 10     i=1,n
!xx   10 SUM = SUM + X(I)                                                  
sum = sum + x(i)
10 continue
!c      XMEAN = SUM / DFLOAT(N)                                           
xmean = sum / dble(n)
do 20   i=1,n
!xx   20 X(I) = X(I) - XMEAN                                               
x(i) = x(i) - xmean
20 continue
!                                                                       
!          VARIANCE COMPUTATION                                         
!                                                                       
sum = 0.d0
do  30     i=1,n
!xx   30 SUM = SUM + X(I) * X(I)                                           
sum = sum + x(i) * x(i)
30 continue
!c      SUM = SUM / DFLOAT(N)                                             
sum = sum / dble(n)
!                                                                       
!          MEAN AND VARIANCE PRINT OUT                                  
!                                                                       
!c      WRITE( 6,7 )   XMEAN, SUM                                         
!
return
!                                                                       
!cx    1 FORMAT( 16I5 )                                                    
!cx    4 FORMAT ( 1H ,10D13.5 )                                            
!cx    5 FORMAT ( 20A4 )                                                   
!cx    6 FORMAT( //1H ,'** ORIGINAL DATA **' )                             
!cx    7 FORMAT ( 1H ,'MEAN =',D15.8,5X,'VARIANCE =',D15.8 )               
!cx    8 FORMAT( 1H ,'TITLE  --  ',1H",20A4,1H" )                          
!cx    9 FORMAT( 1H ,'N =',I5,5X,'FORMAT =',20A4 )                         
!                                                                       
end
!
!
!c       SUBROUTINE  REDUCT( SETX,Z,D,NMK,N0,K,MJ1,LAG,X )                 
subroutine  reduct( setx,z,nmk,n0,k,mj1,lag,x )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          +-----------------------+                                    
!          ! HOUSEHOLDER REDUCTION !                                    
!          +-----------------------+                                    
!                                                                       
!     THIS SUBROUTINE FIRST SETS UP DATA MATRIX X BY AUGMENTING         
!     SUCCESSIVELY SHIFTED ORIGINAL DATA VECTOR Z AND THEN TRANSFORMS   
!     X INTO TRIANGULAR FORM BY HOUSEHOLDER TRANSFORMATION.             
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             HUSHLD                                                    
!             (SETX)                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          SETX:  EXTERNAL SUBROUTINE DESIGNATION                       
!          Z:     ORIGINAL DATA VECTOR                                  
!          D:     WORKING AREA                                          
!          NMK:   DIMENSION OF THE VECTOR OF REGRESSAND (Z(N0+LAG+1),...
!                 Z(N0+LAG+NMK))                                        
!          N0:    INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATION
!          K:     NUMBER OF REGRESSORS                                  
!          MJ1:   ABSOLUTE DIMENSION OF X                               
!          LAG:   MAXIMUM TIME LAG OF THE MODEL                         
!                                                                       
!       OUTPUT:                                                         
!          X:     REDUCED MATRIX ( UPPER TRIANGULAR FORM )              
!                                                                       
!xx      IMPLICIT  REAL*8( A-H,O-Z )                                       
!C      DIMENSION  X(MJ1,1) , D(1)                                        
!C      REAL * 4   Z(1)                                                   
!x      DIMENSION  X(MJ1,1) , Z(1)
!xx      DIMENSION  X(MJ1,LAG+1) , Z(N0+LAG+NMK)
integer nmk, n0, k, mj1, lag
real(dp) z(n0+lag+nmk), x(mj1,lag+1)
! local
integer k1, l, lk, n1, n2
!                                                                       
l = min0( nmk,mj1 )
k1 = k + 1
n1 = l
!
!          +-----------------+                                        +-
!          ! MATRIX X SET UP !                                        ! 
!          +-----------------+                                        +-
!                                                                       
call  setx( z,n0,l,k,mj1,0,lag,x )
!                                                                       
!          +----------------------------+                             +-
!          ! HOUSEHOLDER TRANSFORMATION !                             ! 
!          +----------------------------+                             +-
!                                                                       
!c      CALL  HUSHLD( X,D,MJ1,L,K1 )                                      
call  hushld( x,mj1,l,k1 )
!                                                                       
if( n1 .ge. nmk )     return
!                                                           +------->>--
!                                                           !           
100 l = min0( nmk-n1,mj1-k1 )
!                                                           !           
lk = l + k1
n2 = n0 + n1
!                                                           !           
!          +-----------------+                              !         +-
!          ! MATRIX X SET UP !                              !         ! 
!          +-----------------+                              !         +-
!                                                           !           
call  setx( z,n2,l,k,mj1,1,lag,x )
!                                                           !           
!          +----------------------------+                   !         +-
!          ! HOUSEHOLDER TRANSFORMATION !                   !         ! 
!          +----------------------------+                   !         +-
!                                                           !           
!c      CALL  HUSHLD( X,D,MJ1,LK,K1 )                                     
call  hushld( x,mj1,lk,k1 )
!                                                           !         +-
!                                                           +<-- NO --!N
!                                                                     +-
n1 = n1 + l
if( n1 .lt. nmk )     go to 100
!                                                                       
!                                                                     +-
!                                                                     ! 
!                                                                     +-
return
!                                                                       
end
!
!
!c      SUBROUTINE  SDCOMP( X,A,Y,N,K,MJ,SD )                             
subroutine  sdcomp( x,a,n,k,mj,sd )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE COMPUTES THE RESIDUAL VARIANCE OF THE REGRESSION M
!     WITH THE REGRESSION COEFFICIENTS A(I) (I=1,K).                    
!                                                                       
!       INPUTS:                                                         
!         X:   N*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE REDUCT   
!         A:   VECTOR OF REGRESSION COEFFICIENTS                        
!         Y:   WORKING AREA                                             
!         N:   DATA LENGTH                                              
!         K:   HIGHEST ORDER OF THE MODEL                               
!         MJ:  ABSOLUTE DIMENSION OF X                                  
!                                                                       
!       OUTPUT:                                                         
!         SD:  RESIDUAL VARIANCE                                        
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  X(MJ,1) , Y(1) , A(1)                                  
!x      DIMENSION  X(MJ,1) , Y(K+1) , A(K)
!xx      DIMENSION  X(MJ,K+1) , Y(K+1) , A(K)
integer n, k, mj
real(dp) x(mj,k+1), a(k), sd
! local
integer i, j, k1
real(dp) y(k+1), sum
!                                                                       
k1 = k + 1
!                                                                       
do  20     i=1,k
sum = 0.d0
do  10     j=i,k
!xx   10 SUM = SUM + X(I,J)*A(J)                                           
sum = sum + x(i,j)*a(j)
10 continue
!xx   20 Y(I) = SUM
y(i) = sum
20 continue
y(k1) = 0.d0
!                                                                       
sum = 0.d0
do  30     i=1,k1
!xx   30 SUM = SUM + (Y(I)-X(I,K1))**2                                     
sum = sum + (y(i)-x(i,k1))**2
30 continue
sd = sum / n
!                                                                       
return
!                                                                       
end
!
!
subroutine  setx1( z,n0,l,k,mj1,jsw,lag,x )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          +-----------------+                                          
!          ! MATRIX X SET UP !                                          
!          +-----------------+                                          
!                                                                       
!     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA VECTOR Z(I) (I=N0
!     N0+LAG+L) FOR AUTOREGRESSIVE MODEL FITTING.  X IS THEN USED AS INP
!     SUBROUTINE HUSHLD.                                                
!                                                                       
!       INPUTS:                                                         
!          Z:     ORIGINAL DATA VECTOR                                  
!          N0:    INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATION
!                 (NEW OBSERVATION STARTS AT N0+LAG+1 AND ENDS AT N0+LAG
!          L:     DIMENSION OF THE VECTOR OF NEW OBSERVATIONS           
!          K:     NUMBER OF REGRESSORS (=LAG OR LAG+1)                  
!          MJ1:   ABSOLUTE DIMENSION OF X                               
!          JSW:   =0   TO CONSTRUCT INITIAL L*(K+1) DATA MATRIX         
!                 =1   TO AUGUMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
!                      L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS   
!          LAG:   MAXIMUM TIME LAG OF THE MODEL                         
!                 =K   CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR   
!                 <K   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSO
!                                                                       
!       OUTPUT:                                                         
!          X:      L*(K+1) MATRIX          IF  JSW = 0                  
!                 (K+1+L)*(K+1) MATRIX     IF  JSW = 1                  
!                                                                       
!C      REAL * 8  X(MJ1,1)                                                
!C      DIMENSION  Z(1)                                                   
!x      REAL * 8  X(MJ1,1) , Z(1)
!xx      REAL * 8  X(MJ1,LAG+1) , Z(N0+LAG+L)
integer n0, l, k, mj1, jsw, lag
real(dp) z(n0+lag+l), x(mj1,lag+1)
! local
integer i, ii, i0, j, jj, jksw, k1, ksw
!                                                                       
ksw = 0
if( k .ne. lag )     ksw = 1
k1 = k + 1
i0 = 0
if( jsw .eq. 1 )     i0 = k1
!xx      DO  10     I=1,L
do  20     i=1,l
ii = i + i0
jj = n0 + lag + i
x(ii,k1) = z(jj)
do  10   j=1,lag
jj = jj - 1
jksw = j + ksw
!xx   10 X(II,JKSW) = Z(JJ)                                                
x(ii,jksw) = z(jj)
10 continue
20 continue
!                                                                       
if( ksw .ne. 1 )     return
!                                                                       
!xx      DO  20     I=1,L                                                  
!xx   20 X(I,1) = 1.D0
x(1:l,1) = 1.d0
return
!                                                                       
end
!
!
subroutine  solve( c,r,id,ii,mj2,mj3,g )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE SOLVES THE MATRIX EQUATION  C*G=R,  WHERE THE MATR
!     IS UPPER TRIANGULAR.                                              
!                                                                       
!       INPUTS:                                                         
!          C:      ID*ID COEFFICIENT MATRIX                             
!          R:      ID*ID MATRIX                                         
!          ID:     NUMBER OF ROWS AND COLUMNS OF C AND R                
!          II:     DESIGNATION OF POSITION WHERE THE SOLUSION TO BE STOR
!          MJ2:    ABSOLUTE DIMENSION OF C AND R                        
!          MJ3:    ABSOLUTE DIMENSION OF G AND H                        
!                                                                       
!       OUTPUT:                                                         
!          G:      SOLUTION OF THE MATRIX EQUATION                      
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  C(MJ2,1) , R(MJ2,1) , G(MJ3,MJ3,1)                     
!xx      DIMENSION  C(MJ2,ID) , R(MJ2,ID) , G(MJ3,MJ3,II)
integer id, ii, mj2, mj3
real(dp) c(mj2,id), r(mj2,id), g(mj3,mj3,ii)
! local
integer i, idm1, ip1, j, jj, l
real(dp) sum
!                                                                       
idm1 = id - 1
do  10   j=1,id
!xx   10 G(J,ID,II) = R(ID,J)/C(ID,ID)
g(j,id,ii) = r(id,j)/c(id,id)
10 continue
!                                                                       
!xx      DO  30   JJ=1,IDM1                                                
do  40   jj=1,idm1
i = id - jj
ip1 = i + 1
do  30   j=1,id
sum = 0.d0
do  20   l=ip1,id
!xx   20 SUM = SUM + G(J,L,II)*C(I,L)
sum = sum + g(j,l,ii)*c(i,l)
20 continue
!xx   30 G(J,I,II) = (R(I,J)-SUM) / C(I,I)                                 
g(j,i,ii) = (r(i,j)-sum) / c(i,i)
30 continue
40 continue
!                                                                       
return
end
!
!
!c      SUBROUTINE  SRCOEF( X,M,K,N,MJ,JND,IPR,A,SD )                     
!x      SUBROUTINE  SRCOEF( X,M,K,N,MJ,JND,IPR,A,SD,AIC,IFG,LU )          
!xx      SUBROUTINE  SRCOEF( X,M,K,N,MJ,JND,IPR,A,SD,AIC ) 
subroutine  srcoef( x,m,k,n,mj,jnd,a,sd,aic )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     SUBSET REGRESSION COEFFICIENTS AND RESIDUAL VARIANCE COMPUTATION. 
!                                                                       
!                                                                       
!       INPUTS:                                                         
!         X:     TRIANGULAR MATRIX                                      
!         M:     NUMBER OF REGRESSORS                                   
!         K:     HEIGHEST ORDER OF THE MODELS                           
!         N:     DATA LENGTH                                            
!        JND(I):   (I=1,...,M)  SPECIFICATION OF I-TH REGRESSOR         
!       OUTPUTS:                                                        
!         A:     REGRESSION COEFFICIENTS                                
!         SD:    INNOVATION VARIANCE                                    
!                                                                       
!xx      IMPLICIT  REAL * 8(A-H,O-Z)                                       
!x      DIMENSION  X(MJ,1) , A(1) , JND(1)                                
!xx      DIMENSION  X(MJ,K+1) , A(M) , JND(M)
integer m, k, n, mj, jnd(m)
real(dp) x(mj,k+1), a(m), sd, aic
! local
integer i, i1, ii, j, k1, l, m1, mm1
real(dp) sum, osd
!
k1 = k + 1
m1 = m + 1
!                                                                       
!                                                                       
!          REGRESSION COEFFICIENTS COMPUTATION                          
!                                                                       
l = jnd(m)
a(m) = x(m,k1) / x(m,l)
mm1 = m - 1
if( mm1 .eq. 0 )     go to  60
do  10     ii=1,mm1
i = m - ii
sum = x(i,k1)
i1 = i + 1
do  20     j=i1,m
l = jnd(j)
!xx   20 SUM = SUM - A(J) * X(I,L)                                         
sum = sum - a(j) * x(i,l)
20 continue
l = jnd(i)
!xx   10 A(I) = SUM / X(I,L)                                               
a(i) = sum / x(i,l)
10 continue
!                                                                       
!                                                                       
!          RESIDUAL VARIANCE AND AIC COMPUTATION                        
!                                                                       
60 continue
sd = 0.0d00
do  30     i=m1,k1
!xx   30 SD = SD + X(I,K1) * X(I,K1)                                       
sd = sd + x(i,k1) * x(i,k1)
30 continue
osd = sd / n
aic = n * dlog( osd ) + 2.0d00 * m
!                                                                       
!                                                                       
!          REGRESSION COEFFICIENTS AND RESIDUAL VARIANCE PRINT OUT      
!                                                                       
!x      IF( IPR .LT. 2 )  RETURN                                          
!x      IF( IFG .EQ. 0 )  RETURN
!c      WRITE( 6,5 )                                                      
!c      WRITE( 6,6 )                                                      
!x      IF( IFG.NE.0 )  WRITE( LU,5 )
!x      IF( IFG.NE.0 )  WRITE( LU,6 )
!x      DO  40     I=1,M                                                  
!x      L = JND(I)                                                        
!c   40 WRITE( 6,7 )     L , A(I)                                         
!c      WRITE( 6,8 )     OSD , M , AIC                                    
!x   40 IF( IFG.NE.0 )  WRITE( LU,7 )     L , A(I)
!x      IF( IFG.NE.0 )  WRITE( LU,8 )     OSD , M , AIC
!                                                                       
return
!c    5 FORMAT( 1H0,10X,'SUBSET REGRESSION COEFFICIENTS' )                
!cx    5 FORMAT( /,11X,'SUBSET REGRESSION COEFFICIENTS' )                
!c    6 FORMAT( 1H ,14X,1HI,12X,'A(I)' )                                  
!cx    6 FORMAT( /,15X,1HI,12X,'A(I)' )                                  
!cx    7 FORMAT( 1H ,10X,I5,5F20.10 )                                      
!c    8 FORMAT( 1H0,10X,'SD  = RESIDUAL VARIANCE    =',D19.12,/,11X,'M   =
!cx    8 FORMAT( /,11X,'SD  = RESIDUAL VARIANCE    =',D19.12,/,11X,'M   =
!cx     1 NUMBER OF PARAMETERS =',I3,/,11X,'AIC =  N*LOG(SD) + 2*M     =', 
!xx     2F15.3 )                                                           
end
!
!
subroutine  triinv( x,m,mj,mj1,y )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!       LOWER TRIANGULAR MATRIX INVERSION                               
!                                                                       
!       INPUTS:                                                         
!          X:    TRIANGULAR MATRIX, DIAGONAL ELEMENTS ARE ASSUMED TO BE 
!          M:    DIMENSION OF MATRIX X                                  
!          MJ:   ABSOLUTE DIMENSION OF X                                
!          MJ1:  ABSOLUTE DIMENSION OF Y                                
!       OUTPUT:                                                         
!          Y:    INVERSE OF X                                           
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  X(MJ,1) , Y(MJ1,1)                                     
!xx      DIMENSION  X(MJ,M) , Y(MJ1,M)
integer m, mj, mj1
real(dp) x(mj,m), y(mj1,m)
! local
integer i, ii, ij, j, j1, jj, mm1
real(dp) sum
! local
!
mm1 = m - 1
!xx      DO  10     I=1,MM1                                                
!xx      DO  10     J=I,M                                                  
!xx   10 Y(I,J) = 0.D0
y(1:mm1,1:m) = 0.d0
do  20     i=1,m
!xx   20 Y(I,I) = 1.D0
y(i,i) = 1.d0
20 continue
!xx      DO  40     J=1,MM1
do  41     j=1,mm1
j1 = j + 1
do  40     i=j1,m
sum = 0.d0
ij = i - j
do  30     ii=1,ij
jj = ii + j - 1
!xx   30 SUM = SUM + X(I,JJ) * Y(JJ,J)
sum = sum + x(i,jj) * y(jj,j)
30 continue
!xx   40 Y(I,J) = -SUM
y(i,j) = -sum
40 continue
41 continue
return
end
!
!c      REAL FUNCTION  DMIN*8( X,N )                                      
!xx      REAL*8 FUNCTION  DMIN( X,N )
real(dp) function dmin( x,n )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!       THIS FUNCTION RETURNS THE MINIMUM VALUE AMONG X(I) (I=1,N).     
!                                                                       
!       INPUTS:                                                         
!          X:     VECTOR                                                
!          N:     DIMENSION OF VECTOR X                                 
!                                                                       
!       OUTPUT:                                                         
!          DMIN:  MINIMUM OF X(I) (I=1,N)                               
!                                                                       
!xx      REAL * 8  X(1)
integer n
real(dp) x(n)
! local
integer i
!                                                                       
dmin = x(1)
do  10   i=2,n
if( dmin .gt. x(i) )         dmin = x(i)
10 continue
return
end

