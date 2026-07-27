! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine mulspef(n,k,lagh1,lagh3,cv,p1,p2,ps,pch1,pch2)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.2.2   MULTIPLE SPECTRUM
!-----------------------------------------------------------------------
!      SUBROUTINE MULSPEF(N,K,LAGH1,LAGH3,IR0,IR1,IR2,IC0,IC1,IC2,
!      SUBROUTINE CROSSP(FC,FS,P1,P2,LAGH1,A,LA1)
!      SUBROUTINE ECORSI(FS,LAGH1,FS1,LAGSHF,LA1)
!      SUBROUTINE FGERSI(G,LGP1,FS,LF1)
!      SUBROUTINE SIMCOH(P1,P2,C,S,P3,LAGH1)
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "DAINAMIKKU SISTEMU NO TOKEI-TEKI KAISEKI TO SEIGYO (STATISTICA
!         ANALYSIS AND CONTROL OF DYNAMIC SYSTEMS)" BY H. AKAIKE AND
!         T. NAKAGAWA, SAIENSU-SHA, TOKYO, 1972 (IN JAPANESE)
!-----------------------------------------------------------------------
!     THIS PROGRAM COMPUTES MULTIPLE SPECTRUM ESTIMATES FROM THE OUTPUT
!     OF PROGRAM 5.1.2 MULCOR, USING WINDOWS W1 AND W2.
!     ONLY ONE CARD OF LAGH(MAXIMUM LAG OF COVARIANCES TO BE USED FOR
!     SPECTRUM COMPUTATION) SHOULD BE ADDED ON TOP OF THE OUTPUT OF
!     PROGRAM 5.1.2 MULCOR TO FORM THE INPUT TO THIS PROGRAM.
!     IN THE CARD OUTPUT OF SPECTRUM MATRIX ON AND LOWER DIAGONAL ARE
!     REAL PARTS AND UPPER DIAGONAL ARE IMAGINARY PARTS OF ON AND LOWER
!     DIAGONAL SPECTRAL ELEMENTS.
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!      DIMENSION C(501),S(501),G(501)                                    
!      DIMENSION FC(501),FS(501),P1(501),P2(501),P3(501)                 
!      DIMENSION A1(10),A2(10)                                           
!      DIMENSION P(10020)                                                
!      DIMENSION X(101,10,10)                                            
!xx      DIMENSION C(LAGH1),S(LAGH1),G(LAGH1),CV(LAGH3,K,K)
!xx      DIMENSION FC(LAGH1),FS(LAGH1)
!xx      DIMENSION P1(LAGH1,K,K),P2(LAGH1,K,K),PS(LAGH1,K)
!xx      DIMENSION PCH1(LAGH1,K,K),PCH2(LAGH1,K,K)
!xx      DIMENSION A1(2),A2(3)
!xx      DIMENSION P(LAGH1*2*K)
integer n, k, lagh1, lagh3
real(dp) cv(lagh3,k,k), p1(lagh1,k,k), p2(lagh1,k,k),&
&ps(lagh1,k), pch1(lagh1,k,k), pch2(lagh1,k,k)
! local
integer i, i1, i2, ii, ii1, im1, isw, iwd, jj, jj1, lag2, mla1,&
&mla2
real(dp) c(lagh1), s(lagh1), g(lagh1), fc(lagh1),&
&fs(lagh1), a1(2), a2(3), p(lagh1*2*k), ct5
!
!     INPUT / OUTPUT DATA FILE OPEN
!      CHARACTER(100) DFNAM
!      DFNAM='mulspe.out'
!      CALL SETWND
!      CALL FLOPN3(DFNAM,NFL)
!      IF (NFL.EQ.0) GO TO 999
!     WINDOW W1 DEFINITION
mla1=2
a1(1)=0.5d-00
a1(2)=0.25d-00
!     WINDOW W2 DEFINITION
mla2=3
a2(1)=0.625d-00
a2(2)=0.25d-00
a2(3)=-0.0625d-00
!     LAGH SPECIFICATION
!      READ(5,1) LAGH
!      LAGH1=LAGH+1
!     READING THE OUTPUTS OF PROGRAM 5.1.2 MULCOR
!      READ(5,1) N,LAGH0,K
!      LAGH3=LAGH0+1
!     INITIAL CONDITION PRINT AND PUNCH OUT
!      WRITE(6,60)
!      WRITE(6,61)
!      WRITE(6,62) N,LAGH,K
!      WRITE(6,63)
!      CALL PRCOL1(A1,1,MLA1,1)
!      WRITE(6,64)
!      CALL PRCOL1(A2,1,MLA2,1)
!      WRITE(7,1) N,LAGH,K
!     COMPUTATION STARTS HERE WITH DISK AS EXTERNAL MEMORY.
!C      REWIND 1
!      IP=0
ct5=0.5d-00
lag2=lagh1+lagh1
do 10 ii=1,k
!     AUTO COVARIANCE INPUT
!      READ(5,1) IR0,IC0
!      READ(5,2) (C(I),I=1,LAGH3)
do 15 i=1,lagh1
!xx   15 C(I)=CV(I,II,II)
c(i)=cv(i,ii,ii)
15 continue
do 20 i=1,lagh1
!xx   20 G(I)=C(I)+C(I)
g(i)=c(i)+c(i)
20 continue
g(1)=ct5*g(1)
g(lagh1)=ct5*g(lagh1)
!     F-COS TRANSFORMATION
call fgerco(g,lagh1,fc,lagh1)
!     SPECTRUM SMOOTHING BY WINDOW W1
!      CALL AUSP(FC,P1,LAGH1,A1,MLA1)
call ausp(fc,p1(1,ii,ii),lagh1,a1,mla1)
!     SPECTRUM SMOOTHING BY WINDOW W2
!      CALL AUSP(FC,P2,LAGH1,A2,MLA2)
call ausp(fc,p2(1,ii,ii),lagh1,a2,mla2)
!     TEST STATISTICS COMPUTATION
!      CALL SIGNIF(P1,P2,P3,LAGH1,N)
call signif(p1(1,ii,ii),p2(1,ii,ii),ps(1,ii),lagh1,n)
!     AUTO SPECTRUM AND TEST STATISTICS PRINT OUT
!      WRITE(6,65) IR0,IC0
!      WRITE(6,65) IR0(II),IC0(II)
!      WRITE(6,66)
!      WRITE(6,67)
!      CALL PRCOL3(P1,P2,P3,1,LAGH1,1)
!      AUTO SPECTRUM STORE
!c	LAG2=LAGH1+LAGH1
ii1=(ii-1)*lag2
do 21 i=1,lagh1
i1=ii1+i
i2=i1+lagh1
!      P(I1)=P1(I)
!   21 P(I2)=P2(I)
p(i1)=p1(i,ii,ii)
!xx   21 P(I2)=P2(I,II,II)
p(i2)=p2(i,ii,ii)
21 continue
if(ii.eq.1) go to 10
!     CROSS COVARIANCE INPUT
im1=ii-1
do 11 jj=1,im1
!      READ(5,1) IR1,IC1
!      READ(5,2) (C(I),I=1,LAGH3)
!      READ(5,1) IR2,IC2
!      READ(5,2) (S(I),I=1,LAGH3)
do 25 i=1,lagh1
c(i)=cv(i,ii,jj)
!xx   25 S(I)=CV(I,JJ,II)
s(i)=cv(i,jj,ii)
25 continue
!     F-COS TRANSFORMATION
do 30 i=1,lagh1
!xx   30 G(I)=C(I)+S(I)
g(i)=c(i)+s(i)
30 continue
g(1)=ct5*g(1)
g(lagh1)=ct5*g(lagh1)
call fgerco(g,lagh1,fc,lagh1)
!     F-SIN TRANSFORMATION
do 31 i=1,lagh1
!xx   31 G(I)=S(I)-C(I)
g(i)=s(i)-c(i)
31 continue
g(1)=ct5*g(1)
g(lagh1)=ct5*g(lagh1)
call fgersi(g,lagh1,fs,lagh1)
!     SMOOTHING BY WINDOW W1
isw=1
iwd=0
!      CALL CROSSP(FC,FS,P1,P2,LAGH1,A1,MLA1)
call crossp(fc,fs,p1(1,ii,jj),p1(1,jj,ii),lagh1,a1,mla1)
!     SIMPLE COHERENCE COMPUTATION
33 ii1=(ii-1)*lag2+iwd
jj1=(jj-1)*lag2+iwd
do 32 i=1,lagh1
i1=ii1+i
i2=jj1+i
c(i)=p(i1)
!xx   32 S(I)=P(I2)
s(i)=p(i2)
32 continue
!      CALL SIMCOH(P1,P2,C,S,P3,LAGH1)
if (isw.eq.1) then
call simcoh(p1(1,ii,jj),p1(1,jj,ii),c,s,pch1(1,ii,jj),lagh1)
else
call simcoh(p2(1,ii,jj),p2(1,jj,ii),c,s,pch2(1,ii,jj),lagh1)
end if
!     CROSS SPECTRUM AND SIMPLE COHERENCE   PRINT OUT
!      WRITE(6,65) IR1,IC1
!      WRITE(6,65) IR1(JJ,II),IC1(JJ,II)
!      IF(ISW.NE.1) GO TO 260
!      WRITE(6,166)
!      GO TO 268
!  260 WRITE(6,266)
!  268 WRITE(6,167)
!      CALL PRCOL3(P1,P2,P3,1,LAGH1,1)
if(isw.lt.0) go to 11
!     CROSS SPECTRUM STORE (DISK)
!C     WRITE(1) (P1(I),I=1,LAGH1),(P2(I),I=1,LAGH1)
!     SMOOTHING BY WINDOW W2
isw=-1
iwd=lagh1
!      CALL CROSSP(FC,FS,P1,P2,LAGH1,A2,MLA2)
call crossp(fc,fs,p2(1,ii,jj),p2(1,jj,ii),lagh1,a2,mla2)
go to 33
11 continue
10 continue
!C      END FILE 1
!     SPECTRUM (SMOOTHED BY WINDOW W1) PUNCH OUT
!      LC=LAGH1
!      IL=101
!      ILM1=IL-1
!      J1=0
!      J2=0
!      IB=0
!  416 J1=J2+1
!      J2=J1+ILM1
!      IF(J2.LE.LC) GO TO 413
!  412 J2=LC
!  413 CONTINUE
!C      REWIND 1
!      IP=0
!      DO 500 II=1,K
!      II1=(II-1)*LAG2
!      DO 520 I=J1,J2
!      I0=I-IB
!      I1=II1+I
!  520 X(I0,II,II)=P(I1)
!      IF(II.EQ.1) GO TO 500
!      IM1=II-1
!      DO 501 JJ=1,IM1
!C      READ(1) (P1(I),I=1,LAGH1),(P2(I),I=1,LAGH1)
!      DO 530 I=J1,J2
!      I0=I-IB
!      X(I0,II,JJ)=P1(I)
!  530 X(I0,JJ,II)=P2(I)
!  501 CONTINUE
!  500 CONTINUE
!      DO 610 I=J1,J2
!      I0=I-IB
!      DO 611 II=1,K
!  611 WRITE(7,2) (X(I0,II,JJ),JJ=1,K)
!  610 CONTINUE
!  417 IB=IB+IL
!      IF(J2.LT.LC) GO TO 416
!      CALL FLCLS3(NFL)
!  999 CONTINUE
!    1 FORMAT(10I5)
!    2 FORMAT(4D20.10)
!   60 FORMAT(1H ,13HPROGRAM 5.2.2,3X,17HMULTIPLE SPECTRUM)
!   61 FORMAT(1H ,17HINITIAL CONDITION)
!   62 FORMAT(1H ,2HN=,I5,5X,5HLAGH=,I5,5X,2HK=,I5)
!   63 FORMAT(1H ,12X,9HWINDOW W1/1H ,4X,1HI,11X,5HA1(I))
!   64 FORMAT(1H ,12X,9HWINDOW W2/1H ,4X,1HI,11X,5HA2(I))
!   65 FORMAT(//1H ,8HP(II,JJ),5X,3HII=,I5,3X,3HJJ=,I5)
!   66 FORMAT(1H ,7X,14HPOWER SPECTRUM)
!   67 FORMAT(1H ,4X,1HI,8X,8HPOWER W1,6X,8HPOWER W2,2X,12HSIGNIFICANCE)
!  166 FORMAT(1H ,7X,14HCROSS SPECTRUM,8X,2HW1)
!  266 FORMAT(1H ,7X,14HCROSS SPECTRUM,8X,2HW2)
!  167 FORMAT(1H ,4X,1HI,5X,11HCO-SPECTRUM,1X,13HQUAD-SPECTRUM,2X,16HSIMP
!     ALE COHERENCE)
return
end subroutine
!
subroutine crossp(fc,fs,p1,p2,lagh1,a,la1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES SMOOTHED CROSS SPECTRUM.
!     FC,FS: OUTPUTS OF FGERCO AND FGERSI
!     P1,P2: REAL AND IMAGINARY PART OF SMOOTHED CROSS SPECTRUM
!     LAGH1: DIMENSION OF FC, FS AND PI (I=1,2)
!     A: SMOOTHING COEFFICIENTS
!     LA1: DIMENSION OF A (LESS THAN 11)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION FC(LAGH1),FS(LAGH1),P1(LAGH1),P2(LAGH1)
!xx      DIMENSION A(LA1)
!xx      DIMENSION FC1(521),FS1(521)
integer lagh1, la1
real(dp) fc(lagh1), fs(lagh1), p1(lagh1), p2(lagh1),&
&a(la1)
! local
integer la, lagshf
real(dp) fc1(521), fs1(521)
la=la1-1
lagshf=lagh1+2*la
!     FC SHIFT-RIGHT BY LA FOR END CORRECTION
call ecorco(fc,lagh1,fc1,lagshf,la1)
!     REAL PART SMOOTHING
call smospe(fc1,lagshf,a,la1,p1,lagh1)
!     FS SHIFT-RIGHT BY LA FOR END CORRECTION
call ecorsi(fs,lagh1,fs1,lagshf,la1)
!     IMAGINARY PART SMOOTHING
call smospe(fs1,lagshf,a,la1,p2,lagh1)
return
end subroutine
!
subroutine ecorsi(fs,lagh1,fs1,lagshf,la1)
  use timsac_kinds, only: dp
  implicit none
!     FS SHIFT-RIGHT BY LA FOR IMAGINARY PART END CORRECTION
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION FS(LAGH1),FS1(LAGSHF)
integer lagh1, lagshf, la1
real(dp) fs(lagh1), fs1(lagshf)
! local
integer i, i1, i2, i3, i4, la, la2, lagh2
lagh2=lagh1+1
la=la1-1
do 100 i=1,lagh1
i1=lagh2-i
i2=i1+la
!xx  100 FS1(I2)=FS(I1)
fs1(i2)=fs(i1)
100 continue
la2=lagh1+la
do 110 i=1,la
i1=la1-i
i2=la1+i
i3=la2-i
i4=la2+i
fs1(i1)=-fs1(i2)
!xx  110 FS1(I4)=-FS1(I3)
fs1(i4)=-fs1(i3)
110 continue
return
end subroutine
!
subroutine fgersi(g,lgp1,fs,lf1)
  use timsac_kinds, only: dp
  implicit none
!     FOURIER TRANSFORM (GOERTZEL METHOD)
!     THIS SUBROUTINE COMPUTES FOURIER TRANSFORM OF G(I),I=0,1,...,LG AT
!     FREQUENCIES K/(2*LF),K=0,1,...,LF AND RETURNS SIN TRANSFORM IN
!     FS(K).
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION G(LGP1),FS(LF1)
integer lgp1 ,lf1
real(dp) g(lgp1), fs(lf1)
! local
integer i, i2, k, lf, lg, lg3, lg4
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
110 pi=3.1415926536
alf=lf
t=pi/alf
do 10 k=1,lf1
ak=k-1
tk=t*ak
ck=dcos(tk)
sk=dsin(tk)
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
12 fs(k)=sk*um1
10 continue
return
end subroutine
!
subroutine simcoh(p1,p2,c,s,p3,lagh1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES SIMPLE COHERENCE.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION P1(LAGH1),P2(LAGH1),C(LAGH1),S(LAGH1),P3(LAGH1)
integer lagh1
real(dp) p1(lagh1), p2(lagh1), c(lagh1), s(lagh1),&
&p3(lagh1)
! local
integer i
do 10 i=1,lagh1
!xx   10 P3(I)=(P1(I)**2+P2(I)**2)/(C(I)*S(I))
p3(i)=(p1(i)**2+p2(i)**2)/(c(i)*s(i))
10 continue
return
end subroutine
