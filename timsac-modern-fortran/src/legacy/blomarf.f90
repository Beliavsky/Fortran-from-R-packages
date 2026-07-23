! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine  blomarf( zs,n,id,c,lag,ns0,kmax,zmean,zvari,bw,aic,a,&
&e,aicb,lks,lke,m )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM  BLOMAR                                                   
!.......................................................................
!.....PLANNED BY H.AKAIKE...............................................
!.....DESIGNED BY H.AKAIKE AND G.KITAGAWA...............................
!.....PROGRAMMED BY G.KITAGAWA AND F.TADA...............................
!.....ADDRESS: THE INSTITUTE OF STATISTICAL MATHEMATICS, 4-6-7 MINAMI-AZ
!..............MINATO-KU, TOKYO 106, JAPAN..............................
!.....DATE OF THE LATEST REVISION:  MAR. 6,1979.........................
!.......................................................................
!.....THIS PROGRAM WAS ORIGINALLY PUBLISHED IN "TIMSAC-78", BY H.AKAIKE,
!.....G.KITAGAWA, E.ARAHATA AND F.TADA, COMPUTER SCIENCE MONOGRAPHS, NO.
!.....THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO, 1979.............
!.......................................................................
!     TIMSAC 78.3.4.                                                    
!     _                  __                 _            __             
!     BAYESIAN METHOD OF LOCALLY STATIONARY MULTIVARIATE AR MODEL FITTIN
!                                                                       
!     THIS PROGRAM LOCALLY FITS MULTI-VARIATE AUTOREGRESSIVE MODELS TO  
!     NON-STATIONARY TIME SERIES BY A BAYESIAN PROCEDURE.               
!                                                                       
!       --------------------------------------------------------------- 
!       REFERENCES:                                                     
!          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
!          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
!          30,B,351-363.                                                
!          H.AKAIKE(1978), "A BAYESIAN EXTENSION OF THE MINIMUM AIC     
!          PROCEDURE OF AUTOREGRESSIVE MODEL FITTING.",  RESEARCH MEMO. 
!          NO. 126, THE INSTITUTE OF STATISTICAL MATHEMATICS; TOKYO.    
!       --------------------------------------------------------------- 
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
!             MRDATA                                                    
!             MNONSB                                                    
!       --------------------------------------------------------------- 
!       INPUTS REQUIRED;                                                
!          MT:    INPUT DEVICE FOR ORIGINAL DATA (MT=5: CARD READER).   
!          LAG:   UPPER LIMIT OF THE ORDER OF AR-MODEL, MUST BE LESS THA
!                 OR EQUAL TO 50.                                       
!          NS:    LENGTH OF BASIC LOCAL SPAN.                           
!          KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR    
!                 =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR
!                                                                       
!            -- THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE MRDATA 
!          TITLE: SPECIFICATION OF DATA                                 
!          N:     DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 1000.      
!          ID:    DIMENSION OF DATA,  MUST BE LESS THAN 6               
!                       < ID*(LAG+1)+KSW MUST BE LESS THAN 101 >        
!          IFM:   INPUT FORMAT                                          
!          FORM:  INPUT DATA FORMAT SPECIFICATION STATEMENT.            
!                 -- EXAMPLE --     (8F10.5)                            
!          C(J):  CALIBRATION CONSTANT FOR CHANNEL J (J=1,ID)           
!          Z(I,J): ORIGINAL DATA                                        
!            -----------------------------------------------------------
!
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4  Z                                                       
!c      DIMENSION  Z(1500,10)                                             
!c      DIMENSION  X(200,100) , D(200)                                    
!c      DIMENSION  A(10,10,50) , B(10,10,50)                              
!c      DIMENSION  G(10,10,50) , H(10,10,50) , E(10,10)                   
!xx      DIMENSION  ZS(N,ID), Z(N,ID), C(ID)
!xx      DIMENSION  ZMEAN(ID), ZVARI(ID)
!xx      DIMENSION  A(ID,ID,LAG,KMAX) , B(ID,ID,LAG)
!xx      DIMENSION  G(ID,ID,LAG) , H(ID,ID,LAG) , E(ID,ID,KMAX)
!xx      DIMENSION  BW(KMAX,KMAX), AIC(KMAX,KMAX)
!xx      DIMENSION  AICB(KMAX), LKS(KMAX), LKE(KMAX)
!xx      DIMENSION  F1(LAG*ID,ID,KMAX), F2(LAG*ID,ID,KMAX)
!xxcxx      DIMENSION  X(NS0+LAG+1,(LAG+1)*ID)
!xx      DIMENSION  X((LAG+1)*ID*2,(LAG+1)*ID)
integer n, id, lag, ns0, kmax, lks(kmax), lke(kmax), m
real(dp) zs(n,id), c(id), zmean(id), zvari(id),&
&bw(kmax,kmax), aic(kmax,kmax), a(id,id,lag,kmax),&
&e(id,id,kmax), aicb(kmax)
! local
integer kc, kd, ksw, l, lk, lk1, mf, mj, mj1, mj3, mx, ns
real(dp) z(n,id), b(id,id,lag), g(id,id,lag),&
&h(id,id,lag), f1(lag*id,id,kmax),&
&f2(lag*id,id,kmax), x((lag+1)*id*2,(lag+1)*id)
!
!       PARAMETERS:                                                     
!          MJ:    ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
!          MJ1:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
!          MJ2:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
!          MJ3:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
!                                                                       
!c      CHARACTER(100) IFLNAM,OFLNAM
!c      CALL FLNAM2( IFLNAM,OFLNAM,NFL )
!c      IF (NFL.EQ.0) GO TO 999
!c      IF (NFL.EQ.2) THEN
!c         OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR)
!c      ELSE
!c         CALL SETWND
!c      END IF
!
!c      MJ = 1500                                                         
!c      MJ1 = 200                                                         
!c      MJ3 = 10                                                          
!c      KMAX = 10                                                         
mj = n
!xx      MJ1 = NS0+LAG+1
mj1 = (lag+1)*id*2
mj3 = id
ksw = 0
!
bw(1:kmax,1:kmax) = 0.0d0
aic(1:kmax,1:kmax) = 0.0d0
a(1:id,1:id,1:lag,1:kmax) = 0.0d0
e(1:id,1:id,1:kmax) = 0.0d0
aicb(1:kmax) = 0.0d0
lks(1:kmax) = 0
lke(1:kmax) = 0
f1(1:lag*id,1:id,1:kmax) = 0.0d0
f2(1:lag*id,1:id,1:kmax) = 0.0d0
x(1:mj1,1:(lag+1)*id) = 0.0d0
!
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      READ( 5,1 )     LAG , NS                                          
ns = ns0
!                                                                       
!c      WRITE( 6,2 )                                                      
!c      WRITE( 6,4 )                                                      
!c      WRITE( 6,3 )     LAG , NS , MT                                    
!                                                                       
!c      CALL  MRDATA( MT,MJ,Z,N,ID )                                      
call mrdata( zs,z,n,id,c,zmean,zvari )
!c      CLOSE( MT )
!                                                                       
l = 0
kd = lag * id + ksw
mx = kd * 2
kc = 0
!                                                                       
!                                                                       
m = 0
111 continue
!xx      M = M+1
lk = l + lag
lk1 = lk + 1
if( lk1 .ge. n )     go to 300
m = m+1
if( n-lk1 .le. ns )     ns = n - lk
if( n-lk1-ns .lt. mx )     ns = n - lk
!                                                                       
!c      CALL MNONSB( Z,X,D,G,H,E,KSW,LAG,L,NS,ID,KMAX,MJ,MJ1,MJ3,A,B,AIC )
!xx      CALL MNONSB( Z,G,H,E(1,1,M),KSW,LAG,L,NS,ID,KMAX,KC,MJ,MJ1,MJ3,
!xx     *             BW(1,M),AIC(1,M),A(1,1,1,M),B,AICB(M),F1,F2 )
call mnonsb( z,x,g,h,e(1,1,m),ksw,lag,l,ns,id,kmax,kc,mj,mj1,&
&mj3,bw(1,m),aic(1,m),a(1,1,1,m),b,aicb(m),f1,f2 )
!                                                                       
l = l + ns
!                                                                       
!c      LKE = LK + NS                                                     
lke(m) = lk + ns
mf = lag
!c      WRITE( 6,13 )                                                     
!c      WRITE( 6,16 )                                                     
!c      WRITE( 6,14 )     LK1 , LKE                                       
lks(m) = lk1
!c      DO 10  I=1,MF                                                     
!c      WRITE( 6,16 )                                                     
!c      DO 10  II=1,ID                                                    
!c      WRITE( 6,15 )     (A(II,JJ,I),JJ=1,ID)                            
!c      IF( II .EQ. 1 )     WRITE( 6,17 )   I                             
!c      IF( II .NE. 1 )     WRITE( 6,21 )                                 
!c   10 CONTINUE                                                          
!c      WRITE( 6,16 )                                                     
!c      WRITE( 6,19 )     MF , AIC                                        
!c      WRITE( 6,16 )                                                     
!c      WRITE( 6,12 )                                                     
!c      DO 20  I=1,ID                                                     
!c      WRITE( 6,11 )     (E(I,J),J=1,ID)                                 
!c   20 WRITE( 6,17 )     I                                               
!c      WRITE( 6,16 )                                                     
!c      WRITE( 6,18 )                                                     
!                                                                       
go to 111
300 continue
!c      GO TO 999                                                          
!                                                                       
!c  900 CONTINUE
!c      WRITE(6,600) IVAR,OFLNAM
!c      GO TO 999
!
!c  910 CONTINUE
!c      IF (NFL.EQ.2) CLOSE( 6 )
!c#ifdef __linux__
!cC	reopen #6 as stdout
!c      IF (NFL.EQ.2) OPEN(6, FILE='/dev/fd/1')
!c#endif
!cC /* __linux__ */
!cWRITE(6,610) IVAR,IFLNAM
!cC
!xx  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
!xx  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
!
!c  999 CONTINUE
return
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( 1H ,'PROGRAM TIMSAC 78.3.4',//,'   BAYESIAN METHOD OF LOCA
!xx     *LLY STATIONARY MULTIVARIATE AR MODEL FITTING;',//,'  < BASIC AUTOR
!xx     *EGRESSIVE MODEL >' )                                              
!xx    3 FORMAT( //,1H ,'  UPPER LIMIT OF THE ORDER  K =',I3,/,'   BASIC LO
!xx     1CAL SPAN LENGTH  NS =',I4,/,'   ORIGINAL DATA INPUT DEVICE  MT =',
!xx     2I3 )                                                              
!xx    4 FORMAT( //1H ,10X,'Z(N) = A1*Z(N-1) + A2*Z(N-2) + ... + AK*Z(N-K) 
!xx     1+ W(N)',/,1H ,'  WHERE',/,11X,'K:     ORDER OF THE MODEL',/,11X,  
!xx     2'W(N):  INNOVATION' )                                             
!xx   11 FORMAT( 1H ,18X,5D15.5 )                                          
!xx   12 FORMAT( 1H ,10X,1H.,6X,'INNOVATION VARIANCE MATRIX',53X,1H. )     
!xx   13 FORMAT( 1H ,//11X,35(1H.),'  CURRENT MODEL  ',35(1H.) )           
!xx   14 FORMAT( 1H ,10X,1H.,6X,'M',7X,'AM(I,J)',30X,'DATA  Z(K,.); K=',I5,
!xx     11H,,I5,7X,1H. )                                                   
!xx   15 FORMAT( 1H ,18X,5F15.8 )                                          
!xx   16 FORMAT( 1H ,10X,1H.,85X,1H. )                                     
!xx   17 FORMAT( 1H+,10X,1H.,I7,78X,1H. )                                  
!xx   18 FORMAT( 1H ,10X,87(1H.) )                                         
!xx   19 FORMAT( 1H ,10X,1H.,6X,'ORDER =',I5,67X,1H.,/,11X,1H.,6X,'AIC =', 
!xx     1 F15.3,59X,1H. )                                                  
!xx   21 FORMAT( 1H+,10X,1H.,85X,1H. )                                     
end
!c      SUBROUTINE  MNONSB( Z,X,D,G,H,E,KSW,LAG,N0,NS,ID,KMAX,MJ,MJ1,MJ3,A
!c     *,B,AICB )                                                         
!xx      SUBROUTINE  MNONSB( Z,G,H,E,KSW,LAG,N0,NS,ID,KMAX1,KC,MJ,MJ1,MJ3,
!xx     *                    C,AIC,A,B,AICB,F1,F2 )
subroutine  mnonsb( z,x,g,h,e,ksw,lag,n0,ns,id,kmax1,kc,mj,mj1,&
&mj3,c,aic,a,b,aicb,f1,f2 )
  use timsac_kinds, only: dp
  implicit none
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             DMIN                                                      
!             BAYSWT                                                    
!             MARCOF                                                    
!             MBYSAR                                                    
!             MREDCT                                                    
!             MSDCOM                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          Z:     ORIGINAL DATA; Z(K,I) (K=1,N) REPRESENTS THE RECORD OF
!                 THE I-TH CHANNEL                                      
!          X:     WORKING AREA                                          
!          D:     WORKING AREA                                          
!          G:     WORKING AREA (PARTIAL AUTOREGRESSION COEFFICIENT MATRI
!                 FORWARD MODEL)                                        
!          H:     WORKING AREA (PARTIAL AUTOREGRESSION COEFFICIENT MATRI
!                 BACKWARD MODEL)                                       
!          E:     WORKING AREA                                          
!          KSW:   =0   CONSTATNT VECTOR IS NOT INCLUDED AS A REGRESSOR  
!                 =1   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSO
!          LAG:   UPPER LIMIT OF THE ORDER OF AR-MODEL                  
!          N0:    INDEX OF THE END POINT OF THE FORMER SPAN             
!          NS:    LENGTH OF BASIC LOCAL SPAN                            
!          ID:    DIMENSION OF DATA                                     
!          KMAX:  MAXIMUM NUMBER OF PRECEDING MODELS STORED             
!          MJ:    ABSOLUTE DIMENSION OF Z IN THE MAIN PROGRAM           
!          MJ1:   ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM           
!          MJ3:   ABSOLUTE DIMENSION OF A AND B IN THE MAIN PROGRAM     
!                                                                       
!       OUTPUTS:                                                        
!          A:     AUTOREGRESSIVE COEFFICIENT MATRIX OF FORWARD MODEL    
!          B:     AUTOREGRESSIVE COEFFICIENT MATRIX OF BACKWARD MODEL   
!          AICB:  AIC OF THE CURRENT MODEL                              
!                                                                       
!xx      IMPLICIT  REAL *8  ( A-H , O-Z )                                  
!C      REAL*4     Z(MJ,1) , F1(100,10,11) , F2(100,10,11)
!c      DIMENSION  X(MJ1,1) , D(1) , A(MJ3,MJ3,1) , B(MJ3,MJ3,1)          
!x      DIMENSION  Z(MJ,1) , F1(LAG*ID,ID,KMAX1) , F2(LAG*ID,ID,KMAX1)
!x      DIMENSION  X(MJ1,(LAG+1)*ID), A(MJ3,MJ3,1) , B(MJ3,MJ3,1)
!x      DIMENSION  G(MJ3,MJ3,1) , H(MJ3,MJ3,1) , E(MJ3,1)                 
!c      DIMENSION  Y(100,10) , AIC(11) , C(11)                            
!xx      DIMENSION  Z(MJ,ID) , F1(LAG*ID,ID,KMAX1) , F2(LAG*ID,ID,KMAX1)
!xx      DIMENSION  X(MJ1,(LAG+1)*ID), A(MJ3,MJ3,LAG) , B(MJ3,MJ3,LAG)
!xx      DIMENSION  G(MJ3,MJ3,LAG) , H(MJ3,MJ3,LAG) , E(MJ3,ID) 
!xx      DIMENSION  AIC(KMAX1) , C(KMAX1)
integer ksw, lag, n0, ns, id, kmax1, kc, mj, mj1, mj3
real(dp) z(mj,id), x(mj1,(lag+1)*id), g(mj3,mj3,lag),&
&h(mj3,mj3,lag), e(mj3,id), c(kmax1), aic(kmax1),&
&a(mj3,mj3,lag), b(mj3,mj3,lag), aicb,&
&f1(lag*id,id,kmax1), f2(lag*id,id,kmax1)
!c      DATA     KC  / 0 /                                                
!
! local
integer i, ii, im, imin, ipr, j, jj, kc1, kd, kmax
!xx      DIMENSION  SD1(LAG+1), AIC1(LAG+1), DIC1(LAG+1)
!xx      DIMENSION  BW1(LAG+1), BW2(LAG)
real(dp) sd1(lag+1), aic1(lag+1), dic1(lag+1),&
&bw1(lag+1), bw2(lag), aicm, ek, sd, sdmin
!                                                                       
kmax = kmax1-1
!c      MJ5 = 100                                                         
ipr = 0
kd = lag * id
!          -----------------------------------------------              
!          NEW DATA LOADING AND HOUSEHOLDER TRANSFORMATION              
!          -----------------------------------------------              
!c      CALL  MREDCT( Z,D,NS,N0,LAG,ID,MJ,MJ1,KSW,X )        
call  mredct( z,ns,n0,lag,id,mj,mj1,ksw,x )
!                                                                       
!          -------------------------------------                        
!          BAYESIAN MODEL FITTED TO THE NEW SPAN                        
!          -------------------------------------                        
!c      CALL  MBYSAR( X,D,NS,LAG,ID,KSW,IPR,MJ1,MJ3,A,B,G,H,E,AICB,EK )   
!xx      CALL  MBYSAR( X,NS,LAG,ID,KSW,IPR,MJ1,MJ3,SD1,AIC1,DIC1,
call  mbysar( x,ns,lag,id,ksw,mj1,mj3,sd1,aic1,dic1,&
&aicm,sdmin,imin,bw1,bw2,a,b,g,h,e,aicb,ek )
!                                                                       
if( kc .eq. 0 )  go to 20
!          -----------------------------                                
!          "PARCOR'S" SHIFTED AND STORED                                
!          -----------------------------                                
kc1 = kc+1
!xx      DO 10  JJ=1,KC
do 12  jj=1,kc
ii = kc1 - jj
!xx        DO 10  I=1,KD                                                   
do 11  i=1,kd
do 10  j=1,id
f1(i,j,ii+1) = f1(i,j,ii)
!xx   10 F2(I,J,II+1) = F2(I,J,II)
f2(i,j,ii+1) = f2(i,j,ii)
10 continue
11 continue
12 continue
20 im = 0
!xx      DO 30  II=1,LAG                                                   
!xx      DO 30  I=1,ID
do 32  ii=1,lag
do 31  i=1,id
im = im+1
do 30  j=1,id
f1(im,j,1) = g(i,j,ii)
!xx   30 F2(IM,J,1) = H(I,J,II)
f2(im,j,1) = h(i,j,ii)
30 continue
31 continue
32 continue
if( kc .eq. 0 )  go to 100
!          ---------------------------------------------------------    
!          PREDICTION ERROR VARIANCES AND AIC'S OF THE FORMER MODELS    
!          ---------------------------------------------------------    
aic(1) = aicb
do 50  jj=1,kc
im = 0
!xx        DO 40  II=1,LAG                                                 
!xx        DO 40  I=1,ID                                                   
do 42  ii=1,lag
do 41  i=1,id
im = im+1
do 40  j=1,id
g(i,j,ii) = f1(im,j,jj+1)
!xx   40   H(I,J,II) = F2(IM,J,JJ+1)                                       
h(i,j,ii) = f2(im,j,jj+1)
40 continue
41 continue
42 continue
!                                                                       
call  marcof( g,h,id,lag,mj3,a,b )
!c      CALL  MSDCOM( X,A,Y,D,NS,LAG,ID,KSW,IPR,MJ1,MJ3,MJ5,E,SD )        
!xx      CALL  MSDCOM( X,A,NS,LAG,ID,KSW,IPR,MJ1,E,SD )        
call  msdcom( x,a,ns,lag,id,ksw,mj1,e,sd )
!xx   50 AIC(JJ+1) = NS*DLOG( SD ) + ID*(ID+1)
aic(jj+1) = ns*dlog( sd ) + id*(id+1)
50 continue
!          ----------------------------------------                     
!          BAYESIAN WEIGHTS OF THE PRECEDING MODELS                     
!          ----------------------------------------                     
!-----------------------------  06/11/01
!cx      AICM = DMIN( AIC,KC )                                             
aicm = aic(1)
do 55  i=1,kc
!xx   55 IF( AIC(I) .LT. AICM )  AICM = AIC(I)
if( aic(i) .lt. aicm )  aicm = aic(i)
55 continue
!-----------------------------
call  bayswt( aic,aicm,kc,2,c )
!c      WRITE( 6,3 )     C(1) , AIC(1)                                    
!c      DO 60  I=2,KC1                                                    
!c      IM1 = I-1                                                         
!c   60 WRITE( 6,4 )     IM1 , C(I) , AIC(I)                              
!                                                                       
!          ------------------------                                     
!          AVERAGING OF THE MODELS                                      
!          -----------------------                                      
ek = ek*c(1)**2
!xx      DO 70  II=1,LAG                                                   
!xx      DO 70  I=1,ID                                                     
do 72  ii=1,lag
do 71  i=1,id
do 70  j=1,id
!xx   70 B(I,J,II) = A(I,J,II)*C(1)                                        
b(i,j,ii) = a(i,j,ii)*c(1)
70 continue
71 continue
72 continue
!xx      DO 80  I=1,KD                                                     
do 81  i=1,kd
do 80  j=1,id
f1(i,j,1) = f1(i,j,1)*c(1)
!xx   80 F2(I,J,1) = F2(I,J,1)*C(1)                                        
f2(i,j,1) = f2(i,j,1)*c(1)
80 continue
81 continue
!                                                                       
!xx      DO 90  JJ=1,KC                                                    
!xx        DO 90  I=1,KD                                                   
do 92  jj=1,kc
do 91  i=1,kd
do 90  j=1,id
f1(i,j,1) = f1(i,j,1) + f1(i,j,jj+1)*c(jj+1)
!xx   90   F2(I,J,1) = F2(I,J,1) + F2(I,J,JJ+1)*C(JJ+1)
f2(i,j,1) = f2(i,j,1) + f2(i,j,jj+1)*c(jj+1)
90 continue
91 continue
92 continue
!          -----------------------------------------                    
!          PREDICTION ERROR VARIANCE MATRIX COMPUTED                    
!          -----------------------------------------                    
100 im = 0
kc = kc + 1
if( kc .gt. kmax )     kc = kmax
!xx      DO 110  II=1,LAG                                                  
!xx      DO 110  I=1,ID                                                    
do 112  ii=1,lag
do 111  i=1,id
im = im+1
do 110  j=1,id
g(i,j,ii) = f1(im,j,1)
!xx  110   H(I,J,II) = F2(IM,J,1)                                          
h(i,j,ii) = f2(im,j,1)
110 continue
111 continue
112 continue
!                                                                       
call  marcof( g,h,id,lag,mj3,a,b )
!c      CALL  MSDCOM( X,A,Y,D,NS,LAG,ID,KSW,IPR,MJ1,MJ3,MJ5,E,SD )        
!xx      CALL  MSDCOM( X,A,NS,LAG,ID,KSW,IPR,MJ1,E,SD )
call  msdcom( x,a,ns,lag,id,ksw,mj1,e,sd )
aicb = ns*dlog( sd ) + 2.d0*(ek + ksw*id) + id*(id+1)
!                                                                       
return
!xx    3 FORMAT( ///1H ,13X,'AR-MODEL FITTED TO  !  BAYESIAN WEIGHTS  ! AIC
!xx     1 WITH RESPECT TO THE PRESENT DATA',/,10X,83(1H-),/,1H ,11X,'CURREN
!xx     2T BLOCK',9X,'!',F13.5,7X,'!',F21.3 )
!xx    4 FORMAT( 1H ,6X,I5,' PERIOD FORMER BLOCK  !',F13.5,7X,'!',F21.3 )
end
