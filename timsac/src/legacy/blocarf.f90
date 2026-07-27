! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine blocarf( zs,n,lag,ns0,kmax,zmean,sum,aic,c,b,a,sd,np,&
&ne,sxx )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM  BLOCAR                                                 
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
!     TIMSAC 78.3.2.                                                    
!     _                  ___                __                          
!     BAYESIAN METHOD OF LOCALLY STATIONARY AR MODEL FITTING; SCALAR CAS
!                                                                       
!     THIS PROGRAM LOCALLY FITS AUTOREGRESSIVE MODELS TO NON-STATIONARY 
!     SERIES BY A BAYESIAN PROCEDURE.  POWER SPECTRA FOR STATIONARY SPAN
!     ARE GRAPHICALLY PRINTED OUT.  (THIS PROGRAM IS TENTATIVE.)        
!                                                                       
!     INPUTS REQUIRED:                                                  
!             MT:       INPUT DEVICE FOR ORIGINAL DATA (MT=5 : CARD READ
!             LAG:      UPPER LIMIT OF THE ORDER OF AR MODEL, MUST BE LE
!                       OR EQUAL TO 50.                                 
!             NS:       LENGTH OF BASIC LOCAL SPAN                      
!             KSW:      =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESS
!                       =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REG
!                                                                       
!               -- THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDA
!             TITLE:    SPECIFICATION OF DATA                           
!             N:        DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000
!             DFORM:    INPUT DATA SPECIFICATION STATEMENT.             
!                       -- EXAMPLE  --     (8F10.5)                     
!             (Z(I),I=1,N):  ORIGINAL DATA                              
!               --------------------------------------------------------
!                                                                       
!                                                                       
!     AT EACH STAGE OF MODELLING OF LOCAL AR MODEL, A TWO-STEP BAYESIAN 
!     PROCEDURE IS APPLIED                                              
!        1.   AVERAGING OF THE MODELS WITH DIFFERENT ORDERS FITTED TO TH
!             NEWLY OBTAINED DATA (BY SUBROUTINE ARBAYS).               
!        2.   AVERAGING OF THE MODELS FITTED TO THE PRESENT AND PRECEDIN
!             (AS FAR AS KMAX(<21)) SPANS.                              
!     AS TO THE AVERAGING 1, SEE THE COMMENTS OF PROGRAM UNIBAR.        
!                                                                       
!     AIC OF THE MODEL FITTED TO THE NEW SPAN IS DEFINED BY             
!                                                                       
!          AIC = NS * LOG( SD ) + 2 * EK,                               
!     WHERE                                                             
!          NS:     LENGTH OF NEW DATA                                   
!          SD:     INNOVATION VARIANCE                                  
!          EK:     EQUIVALENT NUMBER OF PARAMETERS, DEFINED AS THE SUM O
!                  SQUARES OF THE BAYESIAN WEIGHTS                      
!                                                                       
!     AIC OF THE MODELS FITTED TO THE PRECEDING SPANS ARE DEFINED BY    
!                                                                       
!          AIC(J+1) = NS * LOG( SD(J) ) + 2     (J=1,KC)                
!     WHERE                                                             
!          SD(J):  PREDICTION ERROR VARIANCE BY THE MODEL FITTED TO J   
!                  PERIODS FORMER SPAN.                                 
!       --------------------------------------------------------------- 
!       REFERENCES:                                                     
!          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
!          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
!          30,B,351-363.                                                
!                                                                       
!          H.AKAIKE(1978), "A BAYESIAN EXTENSION OF THE MINIMUM AIC     
!          PROCEDURE OF AUTOREGRESSIVE MODEL FITTING.",  RESEARCH MEMO. 
!          NO. 126, THE INSTITUTE OF STATISTICAL MATHEMATICS; TOKYO.    
!       --------------------------------------------------------------- 
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
!             REDATA                                                    
!             NONSTB                                                    
!             PRINTA                                                    
!             NRASPE                                                    
!       --------------------------------------------------------------- 
!                                                                       
!c      !DEC$ ATTRIBUTES DLLEXPORT :: BLOCARF
!
!C      IMPLICIT REAL * 8 ( A-H , O-Y )                                   
!xx      IMPLICIT REAL * 8 ( A-H , O-Z )
!c      REAL * 4  TITLE(20) , TTL(13)                                     
!c      DIMENSION Z(10000)                                                
!c      DIMENSION  X(200,51) , D(200) , A(50)                             
!xx      DIMENSION  Z(N), ZS(N)                                            
!xx      DIMENSION  X(NS0,LAG+1), A(LAG,KMAX)
!xx      DIMENSION  AIC(KMAX,KMAX), C(KMAX,KMAX)
!xx      DIMENSION  B(LAG,KMAX), SD(KMAX)
!xx      DIMENSION  NP(KMAX), NE(KMAX)
integer n, lag, ns0, kmax, np(kmax), ne(kmax)
real(dp) zs(n), zmean, sum, aic(kmax,kmax), c(kmax,kmax),&
&b(lag,kmax), a(lag,kmax), sd(kmax), sxx(121,kmax)
!c      DATA  TTL / 4H  CU,4HRREN,4HT MO,4HDEL ,4H(AVE,4HRAGE,4H BY ,4HTHE
!c     1 ,4HBAYE,4HSIAN,4H WEI,4HGHTS,4H)    /                            
!
!xx      DIMENSION    SXX(121,KMAX)
!xx      DIMENSION    F(LAG,KMAX)
! local
integer isw, k, k3, kc, ksw, lag1, mj1, n0, nr, ns
real(dp) z(n), x(ns0,lag+1), f(lag,kmax), bb
!
external  setx1
!
!c      CHARACTER(100) IFLNAM,OFLNAM
!c      CALL FLNAM2( IFLNAM,OFLNAM,NFL )
!c      IF (NFL.EQ.0) GO TO 999
!c      IF (NFL.EQ.2) THEN
!c         OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR )
!c      ELSE
!c         CALL SETWND
!c      END IF
!
!
ksw = 0
!c      KMAX = 20                                                         
!c      MJ1 = 200                                                         
mj1 = ns0
isw = 0
!                                                                       
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      READ( 5,1 )     LAG , NS                                          
ns = ns0
!                                                                       
!c      WRITE( 6,3 )                                                      
!c      IF( KSW .EQ. 1 )     WRITE( 6,5 )                                 
!c      IF( KSW .NE. 1 )     WRITE( 6,4 )                                 
!c      WRITE( 6,6 )                                                      
!c      WRITE( 6,2 )   LAG , NS , MT                                      
!                                                                       
!          ---------------------                                        
!          ORIGINAL DATA LOADING                                        
!          ---------------------                                        
!                                                                       
!c      CALL  REDATA( Z,N,MT,TITLE )                                      
!c      SUBROUTINE  REDATA( X,N,MT,TITLE )                                
call  redata( zs,z,n,zmean,sum )
!c      CLOSE( MT )
!                                                                       
!     ---  INITIAL SETTING  ---                                         
lag1 = lag+1
n0 = 0
k = lag + ksw
k3 = k*3
!                                                                       
nr = 0
kc = 0
10 continue
nr = nr+1
!                                                                       
!          -----------------------------------                          
!          LOCALLY STATIONARY AR-MODEL FITTING                          
!          -----------------------------------                          
!                                                                       
!c      CALL  NONSTB( SETX1,Z,X,D,LAG,N0,NS,KMAX,KSW,ISW,TITLE,MJ1,A,SD ) 
!x      CALL  NONSTB( SETX1,Z,X,LAG,N0,NS,KMAX,KSW,ISW,MJ1,KC,F,
call  nonstb( setx1,z,n,x,lag,n0,ns,kmax,ksw,isw,mj1,kc,f,&
&aic(1,nr),c(1,nr),b(1,nr),a(1,nr),sd(nr) )
!                                                                       
!c      NP = N0 + LAG + 1                                                 
!c      NE = NP+NS-1                                                      
!c      CALL  PRINTA( A,SD,K,TTL,13,TITLE,NP,NE )                         
np(nr) = n0 + lag + 1
ne(nr) = np(nr)+ns-1
!                                                                       
!     ---  SPECTRUM DISPLAY  ---                                        
!c      CALL  NRASPE( SD,A,B,K,0,120,TITLE )                              
call  nraspe( sd(nr),a(:,nr),[bb],k,0,120,sxx(:,nr) )
!                                                                       
!                                                                       
!     ---  PREPARATION FOR THE NEXT STEP  ---                           
n0 = n0 + ns
if(n0+ns+lag1.gt.n)      ns = n-n0-lag1
if(n-n0-ns-lag1.lt.k3)   ns = n-n0-lag1
if( n0+lag1 .lt. n )     go to 10
!                                                                       
!                                                                       
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
!c      WRITE(6,610) IVAR,IFLNAM
!c      GO TO 999
!
!xx  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//5X,100A)
!xx  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//5X,100A)
!
return
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( ///1H ,'  FITTING UP TO THE ORDER  K =',I3,'  IS TRIED',/,
!xx     1'   BASIC LOCAL SPAN  NS =',I4,/,'   ORIGINAL DATA INPUT DEVICE  M
!xx     2T =',I3 )                                                         
!xx    3 FORMAT( //' PROGRAM TIMSAC 78.3.2',/,'   BAYESIAN METHOD OF LOCALL
!xx     *Y STATIONARY AR MODEL FITTING;   SCALAR CASE',//,'   < BASIC AUTOR
!xx     *EGRESSIVE MODEL >' )                                              
!xx    4 FORMAT( 1H ,10X,'Z(I) = A(1)*Z(I-1) + A(2)*Z(I-2) + ... + A(M)*Z(I
!xx     1-M) + E(I)' )                                                     
!xx    5 FORMAT( 1H ,10X,'Z(I) = A(1) + A(2)*Z(I-1) + ... + A(M+1)*Z(I-M) +
!xx     1 E(I)' )                                                          
!xx    6 FORMAT( 1H ,2X,'WHERE',/,11X,'M:     ORDER OF THE MODEL',/,11X,'E(
!xx     1I):  GAUSSIAN WHITE NOISE WITH MEAN 0  AND  VARIANCE SD(M).' )    
end
!c      SUBROUTINE  NONSTB( SETX,Z,X,D,LAG,N0,NS,KMAX,KSW,ISW,TITLE,MJ1,A,
!c     1SD )                                                              
!x      SUBROUTINE  NONSTB( SETX,Z,X,LAG,N0,NS,KMAX1,KSW,ISW,MJ1,KC,F,
subroutine  nonstb( setx,z,n,x,lag,n0,ns,kmax1,ksw,isw,mj1,kc,f,&
&aic,c,b,a,sd )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     BAYESIAN TYPE NON-STATIONARY AUTOREGRESSIVE MODEL FITTING PROCEDUR
!                                                                       
!     THIS SUBROUTINE FIRST FITS AN AUTOREGRESSIVE MODEL TO THE NEWLY OB
!     DATA SPAN AND THEN PRODUCES A BAYESIAN MODEL BY AVERAGING THE MODE
!     FITTED TO THE PRESENT AND PRECEDING SPANS.                        
!                                                                       
!     BAYESIAN WEIGHT OF THE MODEL FITTED TO I PERIODS FORMER SPAN IS   
!     DEFINED BY                                                        
!          C(I)  =  CONST * P(I) / (I+1)                                
!     WHERE                                                             
!          CONST  =  NORMALIZING CONSTANT                               
!          P(I)  =  EXP( -0.5*AIC(I) )                                  
!          AIC(I)  =  AIC WITH RESPECT TO THE PRESENT DATA.             
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             DMIN                                                      
!             ARBAYS                                                    
!             ARCOEF                                                    
!             BAYSWT                                                    
!             REDUCT                                                    
!             SDCOMP                                                    
!       ----------------------------------------------------------------
!                                                                       
!     INPUTS:                                                           
!        SETX:  EXTERNAL SUBROUTINE DESIGNATION                         
!        Z:     ORIGINAL DATA, OUTPUT OF SUBROUTINE REDATA              
!        X:     WORKING AREA                                            
!        D:     WORKING AREA                                            
!        LAG:   UPPER LIMIT OF ORDER OF AR-MODEL                        
!        N0:    INDEX OF THE END POINT OF THE FORMER SPAN               
!        NS:    LENGTH OF BASIC LOCAL SPAN                              
!        KMAX:  MAXIMUM NUMBER OF PRECEDING MODELS STORED               
!        KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR      
!               =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR  
!        ISW:   PRINT OUT CONTROL                                       
!        TITLE: TITLE OF DATA                                           
!        MJ1:   ABSOLUTE DIMENSION OF X                                 
!                                                                       
!     OUTPUTS:                                                          
!        A:     AR-COEFFICIENTS OF THE CURRENT MODEL                    
!        SD:    INNOVATION VARIANCE OF THE CURRENT MODEL                
!                                                                       
!xx      IMPLICIT  REAL*8  ( A-H,O-Z )                                     
!C      REAL*4  Z , TITLE                                                 
!c      REAL*4  TITLE
!c      DIMENSION  X(MJ1,1) , D(1) , A(1) , Z(1) , TITLE(1)               
!c      DIMENSION  F(50,20) , AIC(21) , C(50) , B(50)                     
!x      DIMENSION  X(MJ1,1) , D(LAG+KSW+1), A(1) , Z(1)
!xx      DIMENSION  X(MJ1,1) , D(LAG+KSW+1), A(LAG+KSW) , Z(N)
!xx      DIMENSION  F(LAG+KSW,KMAX1) , AIC(KMAX1) , C(KMAX1) , B(LAG+KSW)
!xx      DIMENSION  SDD(LAG+KSW+1), AICC(LAG+KSW+1), DIC(LAG+KSW+1)
!xx      DIMENSION  B1(LAG+KSW), W(LAG+KSW+1)
integer n, lag, n0, ns, kmax1, ksw, isw, mj1, kc
real(dp) z(n), x(mj1,1), f(lag+ksw,kmax1), aic(kmax1),&
&c(kmax1), b(lag+ksw), a(lag+ksw), sd
! local
integer i, ii, imin, j, k, kc1, kmax
real(dp) d(lag+ksw+1), sdd(lag+ksw+1), aicc(lag+ksw+1),&
&dic(lag+ksw+1), b1(lag+ksw), w(lag+ksw+1), aicb,&
&aicm, sdmin, pn
!c      DATA  KC / 0 /                                                    
external  setx
!                                                                       
k = lag + ksw
kmax = kmax1-1
!                                                                       
!     ---  HOUSEHOLDER REDUCTION  ---                                   
!c      CALL  REDUCT( SETX,Z,D,NS,N0,K,MJ1,LAG,X )                        
call  reduct( setx,z,ns,n0,k,mj1,lag,x )
!                                                                       
!                                                                       
!     ---  BAYESIAN MODEL FITTED TO THE NEW SPAN  ---                   
!c      CALL  ARBAYS( X,D,K,LAG,NS,ISW,TITLE,MJ1,A,B,SD,AICB )            
!xx      CALL ARBAYS( X,D,K,LAG,NS,ISW,MJ1,SDD,AICC,DIC,AICM,SDMIN,IMIN,
call arbays( x,d,k,ns,isw,mj1,sdd,aicc,dic,aicm,sdmin,imin,&
&a,b1,b,w,sd,pn,aicb )
!                                                                       
if( kc .eq. 0 )  go to 110
!                                                                       
!                                                                       
!     ---  PREDICTION ERROR VARIANCE AND AIC OF THE FORMER MODELS  ---  
aic(1) = aicb
do 30  j=1,kc
do 20  i=1,k
!xx   20    D(I) = F(I,J)                                                  
d(i) = f(i,j)
20 continue
call  arcoef( d,k,a )
!c        CALL  SDCOMP( X,A,D,NS,K,MJ1,SD )                               
call  sdcomp( x,a,ns,k,mj1,sd )
!xx   30 AIC(J+1) = NS*DLOG( SD ) + 2.D0                                   
aic(j+1) = ns*dlog( sd ) + 2.d0
30 continue
!                                                                       
!                                                                       
!     ---  BAYESIAN WEIGHTS OF THE MODEL  ---                           
!-------------------------------   06/11/01
!cx      AICM = DMIN( AIC,KC )                                             
aicm = aic(1)
do 33  i=1,kc
!xx   33 IF( AIC(I) .LT. AICM )  AICM = AIC(I)
if( aic(i) .lt. aicm )  aicm = aic(i)
33 continue
!-------------------------------
call  bayswt( aic,aicm,kc,2,c )
!                                                                       
kc1 = kc+1
!c      WRITE( 6,3 )     C(1) , AIC(1)                                    
!c      DO 35  I=2,KC1                                                    
!c         IM1 = I-1                                                      
!c   35 WRITE( 6,4 )     IM1 , C(I) , AIC(I)                              
!                                                                       
!                                                                       
!     ---  AVERAGING OF THE MODELS  ---                                 
do 40  i=1,k
!xx   40 B(I) = B(I)*C(1)                                                  
b(i) = b(i)*c(1)
40 continue
do 70  j=1,kc
do 50  i=1,k
!xx   50    A(I) = F(I,J)                                                  
a(i) = f(i,j)
50 continue
do 60  i=1,k
!xx   60    B(I) = B(I) + C(J+1)*A(I) 
b(i) = b(i) + c(j+1)*a(i)
60 continue
70 continue
!c      WRITE( 6,5 )     (B(I),I=1,K)                                     
!                                                                       
!                                                                       
!     ---  AR-COEFFICIENTS OF THE CURRENT BAYESIAN MODEL  ---           
call  arcoef( b,k,a )
!                                                                       
!                                                                       
!     ---  "PARCOR'S" STORED  ---                                       
!xx      DO 100  J=1,KC                                                    
do 101  j=1,kc
ii = kc1-j
do 100  i=1,k
!xx  100 F(I,II+1) = F(I,II)
f(i,ii+1) = f(i,ii)
100 continue
101 continue
110 continue
do 120  i=1,k
!xx  120 F(I,1) = B(I)                                                     
f(i,1) = b(i)
120 continue
kc = min0( kc+1,kmax )
!                                                                       
!                                                                       
!     ---  PREDICTION ERROR VARIANCE COMPUTED  ---                      
!c      CALL  SDCOMP( X,A,D,NS,K,MJ1,SD )                                 
call  sdcomp( x,a,ns,k,mj1,sd )
!                                                                       
return
!                                                                       
!xx    3 FORMAT( ///1H ,13X,'AR-MODEL FITTED TO  !  BAYESIAN WEIGHTS  ! AIC
!xx     1 WITH RESPECT TO THE PRESENT DATA',/,10X,83(1H-),/,1H ,11X,'CURREN
!xx     2T BLOCK',9X,'!',F13.5,7X,'!',F21.3 )                              
!xx    4 FORMAT( 1H ,6X,I5,' PERIOD FORMER BLOCK  !',F13.5,7X,'!',F21.3 )  
!xx    5 FORMAT( //1H ,'PARTIAL AUTOCORRELATION  B(I) (I=1,K)',/,(1X,10D13.
!xx     15) )                                                              
!                                                                       
end
