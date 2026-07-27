! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

!xx      SUBROUTINE BSUBSTF( ZS,N,IMODEL,LAG,K,IL,LG1,LG2,F,CNST,ZMEAN,SUM,
subroutine bsubstf( zs,n,imodel,lag,k,il,lg1,lg2,zmean,sum,&
&m,aicm,sdm,a1,sd,aic,dic,aicb,sdb,ek,a2,ind,c,c1,c2,b,oeic,esum,&
&omean,om,e,emean,vari,skew,peak,cov,sxx )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM BSUBST
!.......................................................................
!.....PLANNED BY H.AKAIKE...............................................
!.....DESIGNED BY H.AKAIKE AND G.KITAGAWA...............................
!.....PROGRAMMED BY G.KITAGAWA AND F.TADA...............................
!.....ADDRESS: THE INSTITUTE OF STATISTICAL MATHEMATICS, 4-6-7 MINAMI-AZ
!..............MINATO-KU, TOKYO 106, JAPAN..............................
!.....DATE OF THE LATEST REVISION:  MAY  14, 1979.......................
!.......................................................................
!.....THIS PROGRAM WAS ORIGINALLY PUBLISHED IN "TIMSAC-78", BY H.AKAIKE,
!.....G.KITAGAWA, E.ARAHATA AND F.TADA, COMPUTER SCIENCE MONOGRAPHS, NO.
!.....THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO, 1979.............
!.......................................................................
!     TIMSAC 78.1.3.                                                    
!     _                 ____ _                                          
!     BAYESIAN TYPE ALL SUBSET ANALYSIS OF TIME SERIES BY A MODEL LINEAR
!     IN PARAMETERS                                                     
!                                                                       
!     THIS PROGRAM PRODUCES BAYESIAN ESTIMATES OF TIME SERIES MODELS SUC
!     PURE AR MODELS, AR-MODELS WITH NON-LINEAR TERMS, AR-MODELS WITH PO
!     TYPE MEAN VALUE FUNCTIONS, ETC.  THE GOODNESS OF FIT OF A MODEL IS
!     CHECKED BY THE ANALYSIS OF SEVERAL STEPS AHEAD PREDICTION ERRORS. 
!     BY PREPARING AN EXTERNAL SUBROUTINE SETX PROPERLY, ANY TIME SERIES
!     WHICH IS LINEAR IN PARAMETERS CAN BE TREATED.                     
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
!             REDATA                                                    
!             REDLAG                                                    
!             SETLAG                                                    
!             REDREG                                                    
!             REDUCT                                                    
!             ARMFIT                                                    
!             SBBAYS                                                    
!             CHECK                                                     
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS REQUIRED:                                                
!          MT:    ORIGINAL DATA INPUT DEVICE SPECIFICATION              
!          IMODEL:=1  AUTOREGRESSIVE MODEL                              
!                 =2  POLYNOMIAL TYPE NON-LINEAR MODEL (LAG'S READ IN ) 
!                 =3  POLYNOMIAL TYPE NON-LINEAR MODEL (LAG'S AUTOMATICA
!                 =4  AR-MODEL WITH POLYNOMIAL MEAN VALUE FUNCTION      
!                 =5  ANY NON-LINEAR MODEL                              
!                 =6  POLYNOMIAL TYPE EXPONENTIALLY DAMPED NON-LINEAR MO
!                 =7  THIS MODEL IS RESERVED FOR THE USER'S OPTIONAL USE
!          LAG:   MAXIMUM TIME LAG USED IN THE MODEL                    
!          K:     NUMBER OF REGRESSORS                                  
!          IL:    PREDICTION ERRORS CHECKING (UP TO IL-STEPS AHEAD) IS R
!                 N*IL SHOULD BE LESS THAN OR EQUAL TO 20000            
!                                                                       
!       --   THE FOLLOWING INPUTS ARE REQUIRED AT SUBROUTINE REDATA   --
!                                                                       
!          TITLE:   ORIGINAL DATA SPECIFICATION                         
!          N:       DATA LENGTH                                         
!          DFORM:   INPUT DATA FORMAT SPECIFICATION STATEMENT           
!                   -- EXAMPLE --  (8F10.5 )                            
!          X(I) (I=1,N):   ORIGINAL DATA                                
!       ----------------------------------------------------------------
!                                                                       
!xx      PARAMETER  ( MJ2 = 101 )
integer, parameter :: mj2 = 101
!
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL*4     Z(10000) , TITLE(20) , TTL(5) , E(20000)               
!c      REAL*4     TITLE(20) , TTL(5)
!c      DIMENSION  Z(10000), E(20000) 
!c      DIMENSION  X(200,101), D(200) , A(100) , B(100)                  
!c      DATA  TTL  / 4H   B,4HAYES,4HIAN ,4HMODE,4HL    /                 
!xx      COMMON     / BBB / L1(50) , L2(50) , L3(50) , SD0 , CONST         
!c      REAL * 4   CONST , SD0                                            
!xx      DIMENSION  ZS(N), Z(N), E(N,IL) 
!xx      DIMENSION  X(N,K+1), A1(K), A2(K), B(K)                  
!xx      DIMENSION  LG1(3,K), LG2(5)
!xx      DIMENSION  SD(K+1), AIC(K+1), DIC(K+1)
!xx      DIMENSION  IND(K), C(K), C1(K+1), C2(K), ESUM(K+1)
!xx      DIMENSION  EMEAN(IL), VARI(IL), SKEW(IL), PEAK(IL), COV(MJ2)
!xx      DIMENSION  SXX(121), FA(N-LAG,IL)
!c      CHARACTER*4 F(20,K+1)
!xx      INTEGER*1  F(80,K+1)
!xx      CHARACTER(4) G
integer n, imodel, lag, k, il, lg1(3,k), lg2(5), m, ind(k)
real(dp) zs(n), zmean, sum, aicm, sdm, a1(k), sd(k+1),&
&aic(k+1), dic(k+1), aicb, sdb, ek, a2(k), c(k),&
&c1(k+1), c2(k), b(k), oeic, esum(k+1), omean,&
&om, e(n,il), emean(il), vari(il), skew(il),&
&peak(il), cov(mj2), sxx(121)
! local
integer i, ipr, isw, kk, lag1, mj, mj1, nmk
real(dp) z(n), x(n,k+1), fa(n-lag,il), sd0, const
character(4) g
!
integer l1, l2, l3, nn, nps
common     / bbb / l1(50) , l2(50) , l3(50) , sd0 , const
common     / eee / g(20,31)
common     / aaa / nn
!                                                                       
!        EXTERNAL SUBROUTINE DECLARATION:                               
!                                                                       
external  setx1
external  setx2
external  setx4
!x      EXTERNAL  SETX5                                                   
!x      EXTERNAL  SETX6                                                   
!c      EXTERNAL  SETX7                                                   
external  prdct1
external  prdct2
!x      EXTERNAL  PRDCT3                                                  
!x      EXTERNAL  PRDCT6                                                  
!
!c	CHARACTER(100) IFLNAM,OFLNAM
!c	CALL FLNAM2( IFLNAM,OFLNAM,NFL )
!c	IF ( NFL.EQ.0 ) GO TO 999
!c	IF ( NFL.EQ.2 ) THEN
!c	   OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR)
!c	ELSE
!c	   CALL SETWND
!c	END IF
!
!
!        PARAMETERS:                                                    
!             MJ1:  ABSOLUTE DIMENSION FOR SUBROUTINE CALL              
!                                                                       
!c      IF ((IMODEL.LE.0) .OR. (IMODEL.GE.8)) GO TO 150
if ((imodel.le.0) .or. (imodel.ge.7)) go to 150
!c      MJ = 1000                                                         
!c      MJ1 = 200                                                         
nn = n
mj = n
mj1 = n
isw = 1
ipr = 2
!                                                                       
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      READ( 5,1 )     IMODEL , LAG , K , IL                             
!c      WRITE( 6,3 )                                                      
!c      IF( IMODEL.EQ.1 )  WRITE( 6,4 )   IMODEL                          
!c      IF( IMODEL.EQ.2 )  WRITE( 6,5 )   IMODEL                          
!c      IF( IMODEL.EQ.3 )  WRITE( 6,5 )   IMODEL                          
!c      IF( IMODEL.EQ.5 )  WRITE( 6,5 )   IMODEL                          
!c      IF( IMODEL.EQ.6 )  WRITE( 6,11 )   IMODEL                         
!c      WRITE( 6,6 )                                                      
!c      IF( IMODEL.EQ.1 )  WRITE( 6,7 )  K                                
!c      IF( IMODEL.EQ.2 )  WRITE( 6,8 )   LAG , K                         
!c      WRITE( 6,2 )     MT                                               
!                                                                       
!                                                                       
!          ---------------------------------------                      
!          ORIGINAL DATA LOADING AND MEAN DELETION                      
!          ---------------------------------------                      
!                                                                       
!c      CALL  REDATA( Z,N,MT,TITLE )                                      
call  redata( zs,z,n,zmean,sum )
nmk = n - lag
lag1 = lag + 1
!                                                                       
!          ---------------------                                        
!          HOUSEHOLDER REDUCTION                                        
!          ---------------------                                        
!                                                                       
!c      GO TO ( 10,20,30,40,50,60,70 ), IMODEL                            
!x      GO TO ( 10,20,30,40,50,60 ), IMODEL                            
!xx      GO TO ( 10,20,30,40 ), IMODEL                            
if ( imodel .eq. 1 ) go to 10
if ( imodel .eq. 2 ) go to 20
if ( imodel .eq. 3 ) go to 30
if ( imodel .eq. 4 ) go to 40
!                                                                       
10 k = lag
!c      CALL  REDUCT( SETX1,Z,D,NMK,0,K,MJ1,LAG,X )                       
call  reduct( setx1,z,nmk,0,k,mj1,lag,x )
go to 100
!                                                                       
!c   20 CALL  REDLAG( K )                                                 
20 continue
do 21 i=1,k
l1(i) = lg1(1,i)
l2(i) = lg1(2,i)
l3(i) = lg1(3,i)
21 continue
!c      CALL  REDUCT( SETX2,Z,D,NMK,0,K,MJ1,LAG,X )                       
call  reduct( setx2,z,nmk,0,k,mj1,lag,x )
go to 100
!                                                                       
!c   30 CALL  SETLAG( K )                                                 
!c      CALL  REDUCT( SETX2,Z,D,NMK,0,K,MJ1,LAG,X )                       
!xx   30 CALL  SETLAG( K,LG2(1),LG2(2),LG2(3),LG2(4),LG2(5) )
30 call  setlag( kk,lg2(1),lg2(2),lg2(3),lg2(4),lg2(5) )
do 31 i=1,k
lg1(1,i) = l1(i)
lg1(2,i) = l2(i)
lg1(3,i) = l3(i)
31 continue
call  reduct( setx2,z,nmk,0,k,mj1,lag,x )
go to 100
!                                                                       
!c   40 CALL  REDUCT( SETX4,Z,D,NMK,0,K,MJ1,LAG,X )                       
40 call  reduct( setx4,z,nmk,0,k,mj1,lag,x )
go to 100
!                                                                       
!c   50 CALL  REDREG( K )                                                 
!x   50 CONTINUE
!x      DO 51 J=1,K+1
!x         DO 51 I=1,20
!x            II = (I-1)*4+1
!x            G(I,J) = CHAR(F(II,J)) // CHAR(F(II+1,J))
!x     *                             // CHAR(F(II+2,J)) //CHAR(F(II+3,J))
!x   51 CONTINUE
!c      CALL  REDUCT( SETX5,Z,D,NMK,0,K,MJ1,LAG,X )                       
!x      CALL  REDUCT( SETX5,Z,NMK,0,K,MJ1,LAG,X )                       
!x      GO TO 100                                                         
!                                                                       
!x   60 CONTINUE                                                          
!x	SD0 = SUM
!x      CONST = CNST
!c      SD0 = 0.D0                                                        
!c      DO 65  I=1,N                                                      
!c   65 SD0 = SD0 + Z(I)*Z(I)                                             
!c      SD0 = SD0 / N                                                     
!c      READ( 5,9 )     CONST                                             
!c      WRITE( 6,12 )   SD0 , CONST                                       
!c      CALL  REDLAG( K )                                                 
!x      DO 66 I=1,K
!x         L1(I) = LG1(1,I)
!x         L2(I) = LG1(2,I)
!x         L3(I) = LG1(3,I)
!x   66 CONTINUE
!c      CALL  REDUCT( SETX6,Z,D,NMK,0,K,MJ1,LAG,X )                       
!x      CALL  REDUCT( SETX6,Z,NMK,0,K,MJ1,LAG,X )                       
!x      GO TO 100                                                         
!                                                                       
!c   70 CALL  REDUCT( SETX7,Z,D,NMK,0,K,MJ1,LAG,X )                       
!c   70 CALL  REDUCT( SETX7,Z,NMK,0,K,MJ1,LAG,X )                       
!                                                                       
100 continue
!c	CLOSE( MT )
!                                                                       
!          ---------------                                              
!          MAICE PROCEDURE                                              
!          ---------------                                              
!                                                                       
!c      CALL  ARMFIT( X,K,LAG,NMK,ISW,TITLE,MJ1,A,SD,M )                  
!x      IFG=0
!x      CALL ARMFIT( X,K,LAG,NMK,ISW,MJ1,A1,M,SD,AIC,DIC,SDM,AICM,
!x     *             IFG,LU )
call armfit( x,k,lag,nmk,isw,mj1,a1,m,sd,aic,dic,sdm,aicm)
!                                                                       
!          ------------------                                           
!          BAYESIAN PROCEDURE                                           
!          ------------------                                           
!                                                                       
!c      CALL  SBBAYS( X,D,K,NMK,IPR,MJ1,A,SD )                            
call  sbbays( x,k,nmk,ipr,mj1,a2,sdb,ek,aicb,ind,c,c1,c2,b,&
&oeic,esum,omean,om  )
if( k.eq.32 ) return
!                                                                       
!c      IF( IMODEL .EQ. 1 )  CALL  NRASPE( SD,A,B,K,0,121,TITLE )         
if( imodel .eq. 1 )  call  nraspe( sdb,a2,b,k,0,120,sxx )
!                                                                       
nps = lag+1
!                                                                       
!          -------------------------                                    
!          PREDICTION ERROR CHECKING                                    
!          -------------------------                                    
!                                                                       
!x      GO TO ( 110,120,120,150,130,140 ), IMODEL                         
!xx      GO TO ( 110,120,120,150 ), IMODEL                         
if ( imodel .eq. 1 ) go to 110
if ( imodel .eq. 2 ) go to 120
if ( imodel .eq. 3 ) go to 120
go to 150
!                                                                       
!c  110 CALL CHECK( PRDCT1,Z,A,K,0,IL,NPS,N,0,MJ,E )                      
!xx  110 CALL CHECK( PRDCT1,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
110 call check( prdct1,z,a2,k,0,il,nps,n,mj,e,fa,emean,vari,skew,&
&peak,cov,mj2 )
go to 150
!                                                                       
!c  120 CALL  CHECK( PRDCT2,Z,A,K,0,IL,NPS,N,0,MJ,E )                     
!xx  120 CALL  CHECK( PRDCT2,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
120 call  check( prdct2,z,a2,k,0,il,nps,n,mj,e,fa,emean,vari,skew,&
&peak,cov,mj2 )
go to 150
!                                                                       
!c  130 CALL CHECK( PRDCT3,Z,A,K,0,IL,NPS,N,0,MJ,E )                      
!x  130 CALL CHECK( PRDCT3,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
!x     *            PEAK,COV,MJ2 )
!                                                                       
!c  140 CALL  CHECK( PRDCT6,Z,A,K,0,IL,NPS,N,0,MJ,E )                     
!x  140 CALL  CHECK( PRDCT6,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
!x     *             PEAK,COV,MJ2 )
!                                                                       
150 continue
return
!                                                                       
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( /1H ,'ORIGINAL DATA INPUT DEVICE  MT =',I4 )              
!xx    3 FORMAT( ' PROGRAM TIMSAC 78.1.3',/,'   SCALAR TIME SERIES MODEL FI
!xx     *TTING;     BAYESIAN PROCEDURE ( ALL SUBSET REGRESSIN TYPE )' )    
!xx    4 FORMAT( //1H ,'MODEL TYPE',I2,/,'   < AUTOREGRESSIVE MODEL >',/,  
!xx     11H ,10X,'Z(I) = A(1)*Z(I-1) + A(2)*Z(I-2) + ... + A(M)*Z(I-M) + E(
!xx     2I)' )                                                             
!xx    5 FORMAT( //1H ,'MODEL TYPE',I2,/,'   < NON-LINEAR MODEL >',/,1H ,10
!xx     1X,'Z(I) = A(1)*Y(I,1) + A(2)*Y(I,2) + ... + A(K)*Y(I,K) + E(I)' ) 
!xx    6 FORMAT( 1H ,2X,'WHERE',/,11X,'M:     ORDER OF THE MODEL',/,11X,'E(
!xx     1I):  GAUSSIAN WHITE NOISE WITH MEAN 0  AND  VARIANCE SD(M).' )    
!xx    7 FORMAT( 1H ,'FITTING UP TO THE ORDER',I3,2X,'IS TRIED' )          
!xx    8 FORMAT( 1H ,'MAXIMUM LAG =',I4,/,' NUMBER OF REGRESSORS =',I4 )   
!xx    9 FORMAT( F10.0 )                                                   
!xx   11 FORMAT( //1H ,'MODEL TYPE',I2,/,'   < EXPONENTIALLY DAMPED NON-LIN
!xx     *EAR MODEL >',/,1H ,10X,'Z(I) = A(1)*Y(I,1) + A(2)*Y(I,2) + ... + A
!xx     *(K)*Y(I,K) + E(I)' )                                              
!xx   12 FORMAT( 1H ,'SIGMA2 =',D15.5,5X,'CONST =',D15.5 )                 
!                                                                       
end
!c      REAL FUNCTION  BICOEF * 8( K,J )                                  
!xx      REAL*8 FUNCTION  BICOEF( K,J )                                  
real(dp) function bicoef( k,j )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS FUNCTION RETURNS BINOMIAL COEFFICIENTS                       
!                                                                       
!          F(K,J) = K]/(J]*(K-J)])                                      
!                                                                       
!       INPUTS:                                                         
!          K:     NUMBER OF OBJECTS                                     
!          J:     NUMBER OF OBJECTS TAKEN                               
!                                                                       
!       OUTPUT:                                                         
!          F:     NUMBER OF COMBINATIONS OF SELECTING J OBJECTS FROM    
!                 A SET OF K OBJECTS                                    
!                                                                       
!xx      IMPLICIT REAL * 8 ( A-H , O-Z )
integer k, j
! local
integer i, kmj
real(dp) sum, di
!                                                                       
kmj = k-j
sum = 0.d0
do 10   i=1,k
di = i
!xx   10 SUM = SUM + DLOG( DI )
sum = sum + dlog( di )
10 continue
!                                                                       
if( j .eq. 0 )   go to 30
do 20   i=1,j
di = i
!xx   20 SUM = SUM - DLOG( DI )
sum = sum - dlog( di )
20 continue
!                                                                       
30 if( kmj .eq. 0 )   go to 50
do 40   i=1,kmj
di = i
!xx   40 SUM = SUM - DLOG( DI )                                            
sum = sum - dlog( di )
40 continue
!                                                                       
50 bicoef = dexp( sum )
return
!                                                                       
end
!c      SUBROUTINE  CHECK( PRDCT,X,A,K,L,IL,NPS,NPE,IPR,MJ,E )            
!xx      SUBROUTINE  CHECK( PRDCT,X,A,K,L,IL,NPS,NPE,IPR,MJ,E,F,EMEAN,VARI,
subroutine  check( prdct,x,a,k,l,il,nps,npe,mj,e,f,emean,vari,&
&skew,peak,cov,mj2 )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE DRAWS HISTGRAMS AND AUTOCOVARIANCE FUNCTION OF ORI
!     DATA OR PREDICTION ERRORS.                                        
!                                                                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             GRAPH1                                                    
!             GRAPH2                                                    
!             MOMENT                                                    
!             (PRDCT)                                                   
!       ----------------------------------------------------------------
!       INPUTS:                                                         
!          PRDCT:  EXTERNAL SUBROUTINE DESIGNATION                      
!          X:      ORIGINAL DATA                                        
!          A:      REGRESSION COEFFICIENTS                              
!          K:      NUMBER OF REGRESSORS                                 
!          L:      MA-ORDER ( THIS ARGUMENT IS ONLY USED FOR THE CHECKIN
!                  OF AR-MA MODEL)                                      
!          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
!                  =0     ANALYSIS OF ORIGINAL DATA                     
!                  >0     ANALYSIS OF MULTI-STEP (UP TO IL) PREDICTION E
!          NPS:    PREDICTION STARTING POSITION                         
!          NPE:    PREDICTION ENDING POSITION                           
!          IPR:    =0  MATRIX OF SEVERAL STEP AHEAD PREDICTION ERRORS SU
!                  =1  MATRIX OF SEVERAL STEP AHEAD PREDICTION ERRORS IS
!                      OUT                                              
!          MJ:     ABSOLUTE DIMENSION OF E IN THE MAIN PROGRAM          
!                                                                       
!       OUTPUT:                                                         
!          E:      SEVERAL-STEPS PREDICTION ERRORS                      
!                                                                       
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4  X , F , E                                               
!x      DIMENSION  X(1) , E(MJ,1)                                         
!c      DIMENSION  A(1) , COV(120)                                        
!xx      DIMENSION  X(NPE) , E(MJ,IL)
!xx      DIMENSION  A(K) , COV(MJ2)                                        
!xx      DIMENSION  EMEAN(IL), VARI(IL), SKEW(IL), PEAK(IL)
!c      COMMON  /COMXX/ F(2000)                                           
!xx      DIMENSION  F(NPE-NPS+1,IL)
integer k, l, il, nps, npe, mj, mj2
real(dp) x(npe), a(k), e(mj,il), f(npe-nps+1,il),&
&emean(il), vari(il), skew(il), peak(il), cov(mj2)
! local
integer i, ie, ii, istep, isw, j, jj, kk, lag1, lagh, n, nmk
real(dp) sum, cov1, sd
!                                                                       
!                                                                       
istep = 1
isw = il
lagh = 100
n = npe - nps - 1
if( lagh .ge. n )     lagh = n - 1
lag1 = lagh + 1
nmk = n - k
if( isw .gt. 0 )     go to 20
!                                                                       
do 10  i=nps,npe
!xx   10 E(I,1) = X(I) 
e(i,1) = x(i)
10 continue
il = 1
go to 36
!                                                                       
!       ---  SEVERAL STEP AHEAD PREDICITON  ---                         
!                                                                       
20 continue
call  prdct( x,a,k,l,il,nps,npe,mj,e )
!                                                                       
!       ---  PREDICTION ERROR  ---                                      
!                                                                       
!xx      DO 30  II=NPS,NPE                                                 
do 31  ii=nps,npe
i = npe-ii+nps
do 30  j=1,il
jj = i-j+1
!xx   30 E(I,J) = X(I) - E(JJ,J)
e(i,j) = x(i) - e(jj,j)
30 continue
31 continue
if( il .eq. 1 )     go to 34
do 35  j=2,il
jj = j-1
!xx      DO 35  I=1,JJ                                                     
do 33  i=1,jj
ii = i+nps-1
!xx   35 E(II,J) = 0.D0                                                    
e(ii,j) = 0.d0
33 continue
35 continue
34 continue
!c      WRITE( 6,690 )                                                    
!c      IF( IPR .EQ. 0 )   GO TO 36                                       
!c      WRITE( 6,670 )     (I,I=1,IL)                                     
!c      DO 37  I=NPS,NPE                                                  
!c   37 WRITE( 6,640 )     I , (E(I,J),J=1,IL)                            
36 continue
!                                                                       
!       ---  MOMENT COMPUTATION  ---                                    
!                                                                       
do 50  kk=1,il
!                                                                       
ii = nps+kk-1
do 40  i=ii,npe
j = i - ii + 1
!c   40 F(J) = E(I,KK)                                                    
!xx   40 F(J,KK) = E(I,KK)
f(j,kk) = e(i,kk)
40 continue
nmk = npe-nps-(kk-2)
!                                                                       
!c      CALL  MOMENT( F,NMK,EMEAN,VARI,SKEW,PEAK )                        
call  moment( f(1,kk),nmk,emean(kk),vari(kk),skew(kk),peak(kk) )
!                                                                       
!c      IF( ISW .GT. 1 )     WRITE( 6,610 )   KK                          
!c      IF( ISW .EQ. 0 )     WRITE( 6,650 )                               
!c      WRITE( 6,620 )     EMEAN , VARI , SKEW , PEAK                     
!                                                                       
!       ---  HISTOGRAM OF F(I)  ---                                     
!                                                                       
!c      SIG = DSQRT(VARI)                                                 
!c      CALL  GRAPH1( F,1,NMK,SIG )                                       
50 continue
!                                                                       
!       ---  AUTOCORRELATION FUNCTION COMPUTATION  ---                  
!                                                                       
do 100  kk=1,il
!                                                                       
do 70   ii=1,lag1
jj = nps + kk - 1
ie = npe - ii + 1
sum = 0.d0
do 60   i=jj,ie
j = i + ii- 1
!xx   60 SUM = SUM + E(I,KK)*E(J,KK)
sum = sum + e(i,kk)*e(j,kk)
60 continue
!xx   70 COV(II) = SUM / (NPE-NPS-KK+2)
cov(ii) = sum / (npe-nps-kk+2)
70 continue
!                                                                       
cov1 = cov(1)
do 80   i=1,lag1
!xx   80 COV(I) = COV(I) / COV1                                            
cov(i) = cov(i) / cov1
80 continue
!                                                                       
!c      IF( ISW .GT. 0 )   WRITE( 6,630 )     KK                          
!c      IF( ISW .EQ. 0 )   WRITE( 6,660 )                                 
!c      WRITE( 6,680 )     (COV(I),I=1,50)                                
!                                                                       
!       ---  AUTOCORRELATION FUNCTION DISPLAY  ---                      
!                                                                       
sd = npe - nps + 1
sd = dsqrt( 1.d0/sd )
!c      CALL  GRAPH2( COV,LAG1,SD )                                       
!                                                                       
if( istep .eq. kk )     go to 110
!                                                                       
100 continue
!                                                                       
!                                                                       
110 continue
!xx      RETURN                                                            
!xx  610 FORMAT( //1H ,I3,'-STEP AHEAD PREDICTION ERROR' )                 
!xx  620 FORMAT( 1H ,'MEAN       =',D15.8,/,' VARIANCE   =',D15.8,/,' SKEWN
!xx     1ESS   =',D15.8,/,' PEAKEDNESS =',D15.8 )                          
!xx  630 FORMAT( //1H ,'AUTOCORRELATION FUNCTION OF ',I3,'-STEP AHEAD PREDI
!xx     1CTION ERROR' )                                                    
!xx  640 FORMAT( 1H ,I5,5X,10D12.4 )                                       
!xx  650 FORMAT( //,' ORIGINAL DATA' )                                     
!xx  660 FORMAT( //1H ,'AUTOCORRELATION FUNCTION OF ORIGINAL DATA' )       
!xx  670 FORMAT( 1H ,10X,10(4X,'J =',I2,3X) )                              
!xx  680 FORMAT( 1H ,10D13.5 )                                             
!xx  690 FORMAT( ///1H ,45(1H-),2X,'<< J-STEP AHEAD PREDICTION ERROR >>',2X
!xx     1,45(1H-) )                                                        
end
subroutine  moment( x,n,f1,f2,f3,f4 )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          +--------------------+                                       
!          ! MOMENT COMPUTATION !                                       
!          +--------------------+                                       
!                                                                       
!     THIS SUBROUTINE COMPUTES MOMENTS.                                 
!                                                                       
!       INPUTS:                                                         
!          X:     ORIGINAL DATA VECTOR                                  
!          N:     DATA LENGTH                                           
!                                                                       
!       OUTPUTS:                                                        
!          F1:    MEAN OF X                                             
!          F2:    VARIANCE OF X                                         
!          F3:    SKEWNESS OF X                                         
!          F4:    PEAKEDNESS OF X                                       
!                                                                       
!xx      IMPLICIT  REAL * 8  ( F )                                         
!C      DIMENSION  X(1)                                                   
!x	REAL*8  X(1)
!xx      REAL*8  X(N)
integer n
real(dp) x(n), f1, f2, f3, f4
! local
integer i
real(dp) fn, fsum, ff
!                                                                       
fn = n
fsum = 0.d0
do  10     i=1,n
!xx   10 FSUM = FSUM + X(I)                                                
fsum = fsum + x(i)
10 continue
!                                                                       
f1 = fsum / fn
!                                                                       
f2 = 0.d0
f3 = 0.d0
f4 = 0.d0
do  20     i=1,n
ff = x(i) - f1
f2 = f2 + ff*ff
f3 = f3 + ff**3
!xx   20 F4 = F4 + FF**4                                                   
f4 = f4 + ff**4
20 continue
!                                                                       
f2 = f2 / fn
f3 = f3 / (fn*f2*dsqrt(f2))
f4 = f4 / (fn*f2*f2)
!                                                                       
return
end
subroutine  prdct1( z,a,m,l,il,nps,npe,mj,ez )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE COMPUTES SEVARAL STEP AHEAD PREDICTION VALUE OF AN
!     AUTOREGRESSIVE MOVING AVERAGE MODEL.                              
!                                                                       
!       INPUTS:                                                         
!          Z:      ORIGINAL DATA VECTOR                                 
!          A:      AR-MA COEFFICIENTS                                   
!          M:      AR-ORDER                                             
!          L:      MA-ORDER                                             
!          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
!          NPS:    PREDICTION STARTING POSITION                         
!          NPE:    PREDICTION ENDING POSITION                           
!          MJ:     ABSOLUTE DIMENSION OF EZ                             
!                                                                       
!       OUTPUT:                                                         
!          EZ:     PREDICTION VALUE MATRIX                              
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4  Z(1) , EZ(MJ,1)                                         
!x      DIMENSION  Z(1) , EZ(MJ,1)
!x	DIMENSION  A(1)
!xx      DIMENSION  Z(NPE) , EZ(MJ,IL)
!xx      DIMENSION  A(M)                                                   
integer m, l, il, nps, npe, mj
real(dp) z(npe), a(m), ez(mj,il)
! local
integer i, i1, i2, ii, ki, kk, kkm1
real(dp) sum
!                                                                       
!                                                                       
do  100     ii=nps,npe
!                                                                       
!xx      DO  90     KK=1,IL                                                
do  91     kk=1,il
kkm1 = kk - 1
sum = 0.d0
if( kk .eq. 1 )     go to 30
do  20     i=1,kkm1
ki = kk - i
!xx   20 SUM = SUM + A(I)*EZ(II,KI)                                        
sum = sum + a(i)*ez(ii,ki)
20 continue
30 if( kk .gt. m )     go to 50
do  40     i=kk,m
i1 = ii + kkm1 - i
!xx   40 SUM = SUM + A(I)*Z(I1)                                            
sum = sum + a(i)*z(i1)
40 continue
!                                                                       
50 if( l .le. 0 )     go to 90
if( kk .gt. l )     go to 90
do  60     i=kk,l
i1 = m + i
i2 = ii + kkm1 - i
if( i2 .ge. ii )     go to 60
sum = sum + a(i1)*(z(i2)-ez(i2,1))
60 continue
90 ez(ii,kk) = sum
91 continue
100 continue
!                                                                       
return
end
subroutine  prdct2( z,a,k,l,il,nps,npe,mj1,ez )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE COMPUTES SEVERAL STEPS AHEAD PREDICTION VALUES OF 
!     NON-LINEAR REGRESSION MODEL.                                      
!                                                                       
!       INPUTS:                                                         
!          Z:      ORIGINAL DATA VECTOR                                 
!          A:      VECTOR OF AR-COEFFICIENTS                            
!          K:      ORDER OF THE AR MODEL                                
!          L:      THIS DUMMY VARIABLE IS NOT REFERENCED IN THIS SUBROUT
!          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
!          NPS:    PREDICTION STARTING POSITION                         
!          NPE:    PREDICTION ENDING POSITION                           
!          MJ1:    ABSOLUTE DIMENSION OF EZ                             
!                                                                       
!       OUTPUT:                                                         
!          EZ:     PREDICTION VALUE MATRIX                              
!                                                                       
!C      IMPLICIT  REAL*8 ( A-D,O-Y )                                      
!xx      IMPLICIT  REAL*8 ( A-H,O-Z ) 
!c      DIMENSION  Z(1) , A(1) , EZ(MJ1,1) , Y(20)                        
!x      DIMENSION  Z(1) , A(1) , EZ(MJ1,1) , Y(IL)                        
!xx      DIMENSION Z(NPE) , A(K) , EZ(MJ1,IL) , Y(IL)
!c      REAL * 4   CSTDMY, SD0DMY
integer k, l, il, nps, npe, mj1
real(dp) z(npe), a(k), ez(mj1,il)
! local
integer i, ii, j, j1, jj,  lag
real(dp) y(il), sum, xx, x
integer lag1, lag2, lag3
real(dp) cstdmy, sd0dmy
common  / bbb /  lag1(50) , lag2(50) , lag3(50), cstdmy, sd0dmy
!C      INTEGER  RETURN                                                   
!                                                                       
!xx      DO 100  II=NPS,NPE                                                
do 110  ii=nps,npe
do 50  j1=1,il
jj = j1-1
sum = 0.d0
do 40  j=1,k
xx = 1.d0
lag = lag1(j)
!C               ASSIGN 10 TO RETURN                                      
!C               GO TO 200                                                
!C   10          XX = XX*X
x = 1.d0
if ( lag .gt. 0 ) then
i = ii+jj-lag
if ( i .ge. ii ) then
i = i - ii + 1
x = y(i)
else
x = z(i)
end if
end if
xx = xx*x
!                                                
lag = lag2(j)
!C               ASSIGN 20 TO RETURN                                      
!C               GO TO 200                                                                                               
!C   20          XX = XX*X
x = 1.d0
if ( lag .gt. 0 ) then
i = ii+jj-lag
if ( i .ge. ii ) then
i = i - ii + 1
x = y(i)
else
x = z(i)
end if
end if
xx = xx*x
!
lag = lag3(j)
!C               ASSIGN 30 TO RETURN                                      
!C               GO TO 200                                                
!C   30          XX = XX*X
x = 1.d0
if ( lag .gt. 0 ) then
i = ii+jj-lag
if ( i .ge. ii ) then
i = i - ii + 1
x = y(i)
else
x = z(i)
end if
end if
xx = xx*x
!                                                                       
!xx   40       SUM = SUM + A(J)*XX 
sum = sum + a(j)*xx
40 continue
!xx   50    Y(J1) = SUM
y(j1) = sum
50 continue
!                                                                       
do  100  j=1,il
ez(ii,j) = y(j)
100 continue
110 continue
!C      GO TO 300
!        -----  INTERNAL SUBROUTINE  -----
!
!C  200 X = 1.D0
!C      IF ( LAG .LE. 0 )     GO TO 220
!C      I = II+JJ-LAG
!C      IF ( I .GE. II )      GO TO 210
!C      X = Z(I)
!C      GO TO 220
!C  210 I = I - II + 1
!C      X = Y(I)
!C  220 GO TO RETURN, ( 10,20,30 )
!        ----------------------------------
!                                                                                                                                                              
!C  300 RETURN
!    L : DUMMY
l = l
return
end
!c      SUBROUTINE  SETLAG( K )                                           
subroutine  setlag( k,lag1,lag2,lag3,lag4,lag5 )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE  PREPARES SPECIFICATION OF REGRESSORS (L1(I),L2(I)
!     (I=1,...,K) FOR THE FITTING OF (POLYNOMIAL TYPE) NON-LINEAR MODEL.
!     THE OUTPUTS ARE USED AS THE INPUTS TO SUBROUTINE SETX2.           
!                                                                       
!       INPUTS:                                                         
!          LAG1:    MAXIMUM TIME LAG OF LINEAR TERM                     
!          LAG2:    MAXIMUM TIME LAG OF SQUARED TERM                    
!          LAG3:    MAXIMUM TIME LAG OF QUADRATIC CROSS TERM            
!          LAG4:    MAXIMUM TIME LAG OF CUBIC TERM                      
!          LAG5:    MAXIMUM TIME LAG OF CUBIC CROSS TERM                
!                                                                       
!       OUTPUTS:                                                        
!          K:       NUMBER OF REGRESSORS                                
!          (L1(I),L2(I),L3(I))  (I=1,K):     SPECIFICATION OF REGRESSORS
!                                                                       
!              ......................................................   
!              I-TH REGRESSOR IS DEFINED BY                             
!                   Z(N-L1(I)) * Z(N-L2(I)) * Z(N-L3(I))                
!              WHERE  0-LAG TERM Z(N-0) IS REPLACED BY THE CONSTANT 1.  
!              ......................................................   
!                                                                       
!c      REAL * 4   CSTDMY,SD0DMY
!xx      REAL * 8   CSTDMY,SD0DMY
integer k, lag1, lag2, lag3, lag4, lag5
! local
integer i, i1, j, l, ll, m
integer l1, l2, l3
real(dp) cstdmy, sd0dmy
common  / bbb /  l1(50), l2(50), l3(50), cstdmy, sd0dmy
!                                                                       
!c      READ( 5,1 )     LAG1,LAG2,LAG3,LAG4,LAG5                          
!c      WRITE( 6,2 )    LAG1,LAG2,LAG3,LAG4,LAG5                          
if(lag1.le.0)  go to 15
do 10  i=1,lag1
l1(i) = i
l2(i) = 0
!xx   10 L3(I) = 0
l3(i) = 0
10 continue
15 k = lag1
!                                                                       
if(lag2.le.0)   go to 30
do 20  i=1,lag2
k = k+1
l1(k) = i
l2(k) = i
!xx   20 L3(K) = 0                                                         
l3(k) = 0
20 continue
!                                                                       
30 if(lag3.le.1)   go to 50
ll = lag3-1
!xx      DO 40  I=1,LL
do 41  i=1,ll
i1 = i+1
do 40  j=i1,lag3
k = k+1
l1(k) = i
l2(k) = j
!xx   40    L3(K) = 0
l3(k) = 0
40 continue
41 continue
50 m  = k
!                                                                       
if(lag4.le.0)   go to 65
do 60  i=1,lag4
k = k+1
l1(k) = i
l2(k) = i
!xx   60 L3(K) = I                                                         
l3(k) = i
60 continue
65 continue
!                                                                       
if(lag5.le.1)   go to 80
!xx      DO 70  I=1,LAG5                                                   
!xx         DO 70  J=I,LAG5                                                
do 72  i=1,lag5
do 71  j=i,lag5
do 70  l=j,lag5
if(i.eq.j .and. j.eq.l)  go to 70
k = k+1
l1(k) = i
l2(k) = j
l3(k) = l
70 continue
71 continue
72 continue
!                                                                       
!c   80 WRITE( 6,3 )                                                      
80 continue
!xx      IF( LAG1 .EQ. 0 )  GO TO 100                                      
!xx      DO 90  I=1,LAG1                                                   
!c   90 WRITE( 6,4 )     I , L1(I)                                        
!xx   90 CONTINUE
!xx  100 J = LAG1+1                                                        
!xx      IF( LAG2+LAG3 .EQ. 0 )   GO TO 120                                
!xx      DO 110  I=J,M                                                     
!c  110 WRITE( 6,5 )     I , L1(I) , L2(I)                                
!xx  110 CONTINUE
!xx  120 J = M+1                                                           
!xx      IF( LAG4+LAG5 .EQ. 0 )   GO TO 140                                
!xx      DO 130  I=J,K                                                     
!c  130 WRITE( 6,6 )     I ,L1(I) , L2(I) ,L3(I)                          
!xx  130 CONTINUE
!xx  140 CONTINUE                                                          
!                                                                       
return
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( /1H ,'LAG1 =',I3,5X,'LAG2 =',I3,5X,'LAG3 =',I3,5X,
!xx     *        'LAG4 =',I3,5X,'LAG5 =',I3 )
!xx    3 FORMAT( 1H ,4X,'M',5X,'REGRESSOR  Y(I,M)' )                       
!xx    4 FORMAT( 1H ,I5,5X,'Z(I-',I2,')' )                                 
!xx    5 FORMAT( 1H ,I5,5X,'Z(I-',I2,') * Z(I-',I2,')' )                   
!xx    6 FORMAT( 1H ,I5,5X,'Z(I-',I2,') * Z(I-',I2,') * Z(I-',I2,')' )     
!                                                                       
end
subroutine  setx2( z,n0,l,k,mj1,jsw,lag,x )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          +----------------------------------------+                   
!          ! MATRIX X SET UP (FOR NON-LINEAR MODEL) !                   
!          +----------------------------------------+                   
!                                                                       
!     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA VECTOR Z(I) (I=N0
!     N0+K+LAG) FOR THE FITTING OF NON-LINEAR AUTOREGRESSIVE MODEL.  X I
!     USED AS THE INPUT TO SUBROUTINE HUSHLD.                           
!                                                                       
!       INPUTS:                                                         
!          Z:      ORIGINAL DATA VECTOR                                 
!          N0:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
!                  (NEW OBSERVATION STARTS AT N0+LAG+1 AND ENDS AT N0+LA
!          L:      DIMENSION OF THE VECTOR OF NEW OBSERVATIONS          
!          K:      NUMBER OF REGRESSORS                                 
!          MJ1:    ABSOLUTE DIMENSION OF X                              
!          JSW:    =0   TO CONSTRUCT INITIAL L*(K+1) DATA MATRIX        
!                  =1   TO AUGMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
!                       L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS  
!          LAG:    MAXIMUM TIME LAG                                     
!          KSW:    THIS DUMMY VARIABLE IS NOT REFERENCED IN THIS SUBROUT
!                                                                       
!--  THE FOLLOWING VARIABLE SPECIFICATION IS GIVEN EITHER BY REDLAG OR S
!         (L1(I) , L2(I) , L3(I))  (I=1,K)                              
!                                                                       
!               I-TH REGRESSOR IS DEFINED BY                            
!                    Z(N-L1(I)) * Z(N-L2(I)) * Z(N-L3(I))               
!               WHERE 0-LAG TERM Z(N-0) IS AUTOMATICALLY REPLACED BY CON
!                                                                       
!       OUTPUT:                                                         
!          X:      L*(K+1) MATRIX           IF  JSW = 0                 
!                  (K+1+L)*(K+1) MATRIX     IF  JSW = 1                 
!                                                                       
!C      REAL * 8  X(MJ1,1)                                                
!C      DIMENSION  Z(1)                                                   
!c      REAL * 4   CSTDMY, SD0DMY
!x      REAL * 8  X(MJ1,1) ,  Z(1) , ZTEM
!xx      REAL * 8  X(MJ1,K+1) ,  Z(N0+LAG+L) , ZTEM
!xx      REAL * 8  CSTDMY, SD0DMY
integer n0, l, k, mj1, jsw, lag
real(dp) z(n0+lag+l), x(mj1,k+1)
! local
integer i, i0, i1, ii, j1, k1, ll1, ll2, ll3, m1, m2, m3
real(dp) ztem
integer l1, l2, l3
real(dp) cstdmy, sd0dmy
common  / bbb /  l1(50) , l2(50) , l3(50), cstdmy, sd0dmy
!                                                                       
k1 = k + 1
i0 = k1*jsw
do  10     i=1,l
i1 = i + i0
j1 = n0 + lag + i
!xx   10 X(I1,K1) = Z(J1)                                                  
x(i1,k1) = z(j1)
10 continue
!                                                                       
do  70     ii=1,k
ll1 = l1(ii)
ll2 = l2(ii)
ll3 = l3(ii)
!xx      DO  60     I=1,L
do  61     i=1,l
ztem = 1.d0
i1 = i + i0
j1 = n0 + lag + i
if( ll1 .eq. 0 )     go to 40
m1 = j1 - ll1
ztem = ztem * z(m1)
40 if( ll2 .eq. 0 )     go to 50
m2 = j1 - ll2
ztem = ztem * z(m2)
50 if( ll3 .eq. 0 )     go to 60
m3 = j1 - ll3
ztem = ztem * z(m3)
60 x(i1,ii) = ztem
61 continue
70 continue
!                                                                       
return
!                                                                       
end
subroutine  setx4( z,no,l,k,mj1,jsw,lag,x )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA VECTOR Z(I) (I=NO
!     NO+K+L) FOR THE FITTING OF AUTOREGRESSIVE MODEL WITH POLYNOMIAL TY
!     VALUE FUNCTION.  X IS THEN USED AS THE INPUT TO SUBROUTINE HUSHLD.
!       INPUTS:                                                         
!          Z:      ORIGINAL DATA VECTOR                                 
!          NO:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
!          L:      DIMENSION OF THE VECTOR OF NEW OBSERVATIONS          
!          K:      NUMBER OF REGRESSORS                                 
!          MJ1:    ABSOLUTE DIMENSION OF X                              
!          JSW:    =0   TO CONSTRUCT INITIAL L*(K+1) DATA MATRIX        
!                  =1   TO AUGMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
!                       L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS  
!          LAG:    MAXIMUM TIME LAG OF THE MODELS                       
!          N:      DATA LENGTH                                          
!                                                                       
!       OUTPUT:                                                         
!          X:      L*(K+1) MATRIX           IF   JSW = 0                
!                  (K+1+L)*(K+1) MATRIX     IF   JSW = 1                
!                                                                       
!                                                                       
!xx      IMPLICIT  REAL  * 8 (A-H,O-Z)                                     
!C      REAL * 4 Z(1)                                                     
!x      DIMENSION  Z(1)
!x      DIMENSION  X(MJ1,1)                                               
!xx      DIMENSION  X(MJ1,K+1) ,  Z(NO+LAG+L)
integer no, l, k, mj1, jsw, lag
real(dp) z(no+lag+l), x(mj1,k+1)
! local
integer i, i0, ii, j, jj, k1, m, m1, nn
real(dp) bn, y, xx
integer n
!c      COMMON     / AAA /  N , M                                         
common  / aaa /  n
!                                                                       
!          M:      ORDER OF POLYNOMIAL OF MEAN VALUE FUNCTION           
!                                                                       
m = k - lag - 1
k1 = k + 1
i0 = jsw*k1
m1 = m + 1
lag=k-m1
bn = 2.d0/(n-lag-1.d0)
!xx      DO 10  I=1,L                                                      
do 11  i=1,l
y= bn*(no+i-1)-1.d0
xx= 1.d0
do 10 j=1,m1
ii= i+i0
x(ii,j) = xx
!xx   10 XX = XX*Y 
xx = xx*y
10 continue
11 continue
!                                                                       
!xx      DO 20  I=1,L 
do 21  i=1,l
ii = i+i0
nn = no+lag+i
x(ii,k1) = z(nn)
do 20   j=1,lag
nn = nn-1
jj = j+m1
!xx   20 X(II,JJ) = Z(NN)                                                  
x(ii,jj) = z(nn)
20 continue
21 continue
!                                                                       
return
!                                                                       
end
subroutine  srtmin( x,n,ix )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!       THIS SUBROUTINE ARRANGES X(I) (I=1,N) IN ORDER OF INCREASING    
!       MAGNITUDE OF X(I)                                               
!                                                                       
!       INPUTS:                                                         
!          X:   VECTOR                                                  
!          N:   DIMENSION OF THE VECTOR                                 
!       OUTPUTS:                                                        
!          X:   ARRANGED VECTOR                                         
!          IND: INDEX OF ARRANGED VECTOR                                
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  X(1) , IX(1)                                           
!xx      DIMENSION  X(N) , IX(N)
integer n, ix(n)
real(dp) x(n)
! local
integer i, ii, it, min, nm1
real(dp) xmin, xt
!                                                                       
nm1 = n - 1
do  30     i=1,n
!xx   30 IX(I) = I                                                         
ix(i) = i
30 continue
do  20     ii=1,nm1
xmin = x(ii)
min = ii
do  10     i=ii,n
if( xmin .lt. x(i) )     go to 10
xmin = x(i)
min = i
10 continue
if( xmin .eq. x(ii) )     go to 20
xt = x(ii)
x(ii) = x(min)
x(min) = xt
it = ix(ii)
ix(ii) = ix(min)
ix(min) = it
20 continue
!                                                                       
return
end
!c      SUBROUTINE  SUBSPC( B,K,N,IPR,EK )                                
!xx      SUBROUTINE  SUBSPC( B,K,N,IPR,EK,IND,C,C1,C2,OEIC,ESUM1,OMEAN,OM )
subroutine  subspc( b,k,n,ek,ind,c,c1,c2,oeic,esum1,omean,om )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!       THIS SUBROUTINE PRODUCES BAYESIAN ESTIMATES OF PARTIAL CORRELATI
!       BY CHECKING ALL SUBSET REGRESSION MODELS.                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             BICOEF                                                    
!             SRTMIN                                                    
!       ----------------------------------------------------------------
!                                                                       
!         INPUTS:                                                       
!           B:   LEAST SQUARES ESTIMATES OF PARTIAL CORRELATIONS        
!           K:   DIMENSION OF VECTOR A                                  
!           N:   NUMBER OF OBSERVATIONS USED FOR THE ESTIMATION OF A(I) 
!           IPR: =0  TO SUPPRESS THE OUTPUTS                            
!                >0  TO PRINT OUT THE OUTPUTS                           
!                                                                       
!         OUTPUTS:                                                      
!           B(I) (I=1,K):   BAYESIAN ESTIMATES OF PARTIAL CORRELATIONS  
!           EK:   EQUIVALENT NUMBER OF FREE PARAMETERS IN THE BAYESIAN M
!                                                                       
!xx      IMPLICIT  REAL * 8 ( A-H , O-Z )                                  
!c      DIMENSION  B(1) , C(50) , D(50,50)                                
!c      DIMENSION  IND(50) , KND(50) , ESUM(50)                           
!c      DIMENSION  C1(50), C2(50), ESUM1(50)
!x      DIMENSION  B(1) , C(K) , D(K+1,K+1)                                
!xx      DIMENSION  B(K) , C(K) , D(K+1,K+1)                                
!xx      DIMENSION  IND(K) , KND(K+1) , ESUM(K+1)                           
!xx      DIMENSION  C1(K+1), C2(K), ESUM1(K+1)
integer k, n, ind(k)
real(dp) b(k), ek, c(k), c1(k+1), c2(k), oeic, esum1(k+1),&
&omean, om
! local
integer i, ip, iq, j, k1, knd(k+1), kmj, m
real(dp) bicoef, d(k+1,k+1), esum(k+1), cc, dn, sum, eic,&
&sumc, exic, osum
!                                                                       
cc = 1.d0 + dlog(2.d0)
k1 = k + 1
dn = n
!xx      DO 10   I=1,K1                                                    
!xx      ESUM(I) = 0.D0                                                    
!xx      DO 10   J=1,K1                                                    
!xx   10 D(I,J) = 0.D0 
esum(1:k1) = 0.d0
d(1:k1,1:k1) = 0.d0
!                                                                       
!          SQUARE OF PARTIAL CORRELATIONS ( NORMALISED BY MULTIPLYING N 
!                                                                       
do 20   i=1,k
!xx   20 C(I) = B(I)*B(I)*DN                                               
c(i) = b(i)*b(i)*dn
20 continue
!                                                                       
!          ARRANGEMENT OF C(I) IN ORDER OF INCREASING MAGNITUDE         
!                                                                       
call  srtmin( c,k,ind )
!c      IF( IPR .LE. 1 )     GO TO 60                                     
!c      WRITE( 6,7 )                                                      
!c      DO  50     I=1,K                                                  
!c   50 WRITE( 6,6 )     I , IND(I) , C(I)                                
!c   60 CONTINUE                                                          
!                                                                       
!          FIND THE MINIMUM OF EIC                                      
!                                                                       
oeic = cc*k
sum = 0.d0
do 30   i=1,k
sum = sum + c(i)
eic = sum + cc*(k-i)
!xx   30 IF( OEIC .GT. EIC )   OEIC = EIC                                  
if( oeic .gt. eic )   oeic = eic
30 continue
!c      WRITE( 6,604 )   OEIC                                             
!                                                                       
!--------  COMPUTATION OF EIC'S OF WHOLE SUBSET REGRESSION MODELS  -----
!                                                                       
!          INITIAL SETTING                                              
!                                                                       
do 40   i=1,k
!xx   40 KND(I) = 0                                                        
knd(i) = 0
40 continue
knd(k1) = 1
sum = 0.d0
sumc = 0.d0
m = k
ip = 0
iq = 0
!                                                                       
100 continue
!                                                                       
!          -----  SPECIFICATION OF NEXT SUBSET  -----                   
!xx               DO 110   I=1,K
do 111   i=1,k
if( knd(i) .eq. 0 )   go to 110
knd(i) = 0
go to 120
110 knd(i) = 1
111 continue
120 continue
!          ------------------------------------------                   
!                                                                       
130 continue
if( ip .gt. k )   go to 200
if( knd(ip+1) .eq. 0 )   go to 140
!                                                                       
if( iq .eq. 0 )   go to 165
if( knd(iq) .eq. 1 )   go to 150
iq = iq-1
sumc = sumc + c(iq+1)
!                                                                       
if( sumc + cc*(k-ip+iq) .gt. oeic + 40.d0 )   go to 180
go to 150
!                                                                       
140 ip = ip+1
iq = ip-1
sumc = c(ip)
if( sumc + cc .gt. oeic + 40.d0 )   go to 200
!                                                                       
150 m = k-ip+iq
sum = sumc
if( iq .eq. 0 )   go to 165
do 160   i=1,iq
if( knd(i) .eq. 1 )   go to 160
m = m-1
sum = sum + c(i)
160 continue
165 continue
eic = sum + cc*m - oeic
if( eic .gt. 40.d0 )   go to 100
exic = dexp( -0.5d0*eic )
esum(m+1) = esum(m+1) + exic
do 170   i=1,k
!xx  170 IF( KND(I) .EQ. 1 )   D(I,M+1) = D(I,M+1) + EXIC                  
if( knd(i) .eq. 1 )   d(i,m+1) = d(i,m+1) + exic
170 continue
go to 100
!         --------------------------------------------                  
180 do 190   i=1,ip
!xx  190          KND(I) = 1                                               
knd(i) = 1
190 continue
knd(ip+1) = 0
ip = ip+1
iq = ip-1
go to 130
!         ---------------------------------------------                 
!                                                                       
!--------------------------  WHOLE SUBSETS CHECKED  --------------------
!                                                                       
200 continue
!c      IF( IPR .GE. 2 )     WRITE( 6,8 )                                 
!c      IF( IPR .GE. 2 )     WRITE( 6,607 )     (ESUM(I),I=1,K1)          
do 201 i=1,k1
!xx  201 ESUM1(I) = ESUM(I)
esum1(i) = esum(i)
201 continue
!                                                                       
!          MEAN OF NUMBER OF PARAMETERS                                 
!                                                                       
osum = 0.d0
sum = esum(1)
do 210   i=1,k
sum = sum + esum(i+1)
!xx  210 OSUM = OSUM + I*ESUM(I+1)                                         
osum = osum + i*esum(i+1)
210 continue
omean = osum / sum
om = omean / k
!c      IF( IPR .GE. 2 )     WRITE( 6,608 )     OMEAN , OM                
!                                                                       
!       --  BINOMIAL TYPE DAMPER  --                                    
!                                                                       
do 220   i=1,k1
j = i-1
kmj = k-j
!c  220 C(I) = BICOEF(K,J)*(OM**J)*((1.D0-OM)**KMJ)                       
!c      C(1) = 1.D0 / (1.D0 + K)                                          
!xx  220 C1(I) = BICOEF(K,J)*(OM**J)*((1.D0-OM)**KMJ)                       
c1(i) = bicoef(k,j)*(om**j)*((1.d0-om)**kmj)
220 continue
c1(1) = 1.d0 / (1.d0 + k)
do 221  i=1,k
!c  221 C(I+1) = C(I) * I / (1.D0 + K - I)                                
!xx  221 C1(I+1) = C1(I) * I / (1.D0 + K - I)                                
c1(i+1) = c1(i) * i / (1.d0 + k - i)
221 continue
!c      IF( IPR .GE. 2 )     WRITE( 6,609 )                               
!c      IF( IPR .GE. 2 )     WRITE( 6,607 )     (C(I),I=1,K1)             
!                                                                       
do 230   i=1,k1
!c  230 ESUM(I) = ESUM(I)*C(I)                                            
!xx  230 ESUM(I) = ESUM(I)*C1(I)                                            
esum(i) = esum(i)*c1(i)
230 continue
!                                                                       
sum = 0.d0
do 240   i=1,k1
!xx  240 SUM = SUM + ESUM(I)                                               
sum = sum + esum(i)
240 continue
!                                                                       
!xx      DO 250   J=1,K1                                                   
do 251   j=1,k1
do 250   i=1,k
!c  250 D(I,J) = D(I,J)*C(J) / SUM                                        
!xx  250 D(I,J) = D(I,J)*C1(J) / SUM                                        
d(i,j) = d(i,j)*c1(j) / sum
250 continue
251 continue
!                                                                       
!          WEIGHTS OF PARTIAL CORRELATIONS                              
!                                                                       
!xx      DO 260   I=1,K                                                    
!c  260 C(I) = 0.D0                                                       
!xx  260 C2(I) = 0.D0
c2(1:k) = 0.d0
!xx      DO 270   I=1,K                                                    
do 271   i=1,k
do 270   j=1,k1
!c  270 C(I) = C(I) + D(I,J)                                              
!xx  270 C2(I) = C2(I) + D(I,J)
c2(i) = c2(i) + d(i,j)
270 continue
271 continue
!c      IF( IPR .GE. 2 )     WRITE( 6,603 )                               
!c      IF( IPR .GE. 2 )     WRITE( 6,607 )     (C(I),I=1,K)              
!                                                                       
!          AVERAGING AND REARRANGEMENT OF PARTIAL CORRELATIONS          
!                                                                       
ek = 1.d0
do 280   i=1,k
j = ind(i)
!c      B(J) = B(J)*C(I)                                                  
!c  280 EK = EK + C(I)**2                                                 
b(j) = b(j)*c2(i)
!xx  280 EK = EK + C2(I)**2                                                 
ek = ek + c2(i)**2
280 continue
!c      IF( IPR .LE. 1 )     RETURN                                       
!c      WRITE( 6,602 )                                                    
!c      DO  290     I=1,K                                                 
!c  290 WRITE( 6,609 )     I , B(I)                                       
!                                                                       
return
!                                                                       
!xx    3 FORMAT( 1H ,20I3 )                                                
!xx    4 FORMAT( 1H ,3D20.10,2I5 )                                         
!xx    6 FORMAT( 1H ,2I7,F13.5 )                                           
!xx    7 FORMAT( 1H ,6X,'I IND(I)',4X,'N*B(I)**2' )                        
!xx    8 FORMAT( 1H ,'ESUM(I) (I=1,M+1)' )                                 
!xx    9 FORMAT( 1H ,'***  BINOMIAL TYPE  ***' )                           
!xx  602 FORMAT( 1H ,'PARTIAL CORRELATIONS OF THE BAYESIAN MODEL',/,1H ,6X,
!xx     1'I',9X,'A(I)' )                                                   
!xx  603 FORMAT( 1H ,'FINAL BAYESIAN WEIGHTS OF PARTIAL CORRELATIONS' )    
!xx  604 FORMAT( 1H ,'  OAIC =',D13.5 )                                    
!xx  606 FORMAT( 1H ,'OSUM =',D20.10,5X,'SUM =',D20.10,5X,'COD =',D20.10 ) 
!xx  607 FORMAT( 1H ,10D13.5 )                                             
!xx  608 FORMAT( 1H ,'OMEAN =',D15.8,5X,'OM =',D15.8 )                     
!xx  609 FORMAT( 1H ,I7,F13.5 )                                            
end
!c      SUBROUTINE  SBBAYS( X,D,K,N,IPR,MJ1,A,SD )                        
subroutine  sbbays( x,k,n,ipr,mj1,a,sd,ek,aic,ind,c,c1,c2,b,&
&oeic,esum,omean,om  )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PRODUCES BAYESIAN MODEL BASED ON ALL SUBSET       
!     REGRESSION MODELS USING THE OUTPUT OF SUBROUTINE REDUCT.          
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             RECOEF                                                    
!             SDCOMP                                                    
!             SUBSPC                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          X:     N*(K+1) TRIANGULAR MATRIX,OUTPUT OF SUBROUTINE REDUCT 
!          K:     NUMBER OF REGRESSORS OF THE BAYESIAN MODEL            
!          N:     DATA LENGTH                                           
!          IPR:   =0  TO SUPPRESS THE OUTPUTS                           
!                 =1  TO PRINT OUT FINAL RESULT                         
!                 =2  TO PRINT OUT INTERIM AND FINAL RESULTS            
!          MJ1:   ABSOLUTE DIMENSION OF X                               
!                                                                       
!       OUTPUTS:                                                        
!          A(I) (I=1,K):   REGRESSION COEFFICIENTS OF BAYESIAN MODEL    
!          SD:    RESIDUAL VARIANCE                                     
!                                                                       
!xx      IMPLICIT  REAL * 8  (A-H , O-Z )                                  
!c      DIMENSION  X(MJ1,1) , A(1) , D(1)                                 
!c      DIMENSION  B(50) , G(50)                                          
!c      DIMENSION  IND(50), C(50), C1(50), C2(50), ESUM(50)
!x      DIMENSION  X(MJ1,1) , A(1) , D(K)                                 
!xx      DIMENSION  X(MJ1,K+1) , A(K) , D(K)                                 
!xx      DIMENSION  B(K) , G(K)
!xx      DIMENSION  IND(K), C(K), C1(K+1), C2(K), ESUM(K+1)
integer k, n, ipr, mj1, ind(k)
real(dp) x(mj1,k+1), a(k), sd, ek, aic, c(k), c1(k+1),&
&c2(k), b(k), oeic, esum(k+1), omean, om
! local
integer i, j, k1
real(dp) d(k), g(k), fn, sum, bb
k1 = k + 1
fn = n
!c      IF( IPR .GE. 2 )     WRITE( 6,3 )                                 
!                                                                       
!          PARTIAL CORRELATIONS COMPUTATION                             
!                                                                       
sum = x(k1,k1)**2
do 10   i=1,k
j = k1-i
sum = sum + x(j,k1)**2
g(j) = dsqrt( sum )
!xx   10 B(J) = X(J,K1)*X(J,J) / (G(J)*DABS(X(J,J)))                       
b(j) = x(j,k1)*x(j,j) / (g(j)*dabs(x(j,j)))
10 continue
!                                                                       
!          PARTIAL CORRELATIONS OF BAYESIAN MODEL COMPUTATION           
!                                                                       
!c      CALL  SUBSPC( B,K,N,IPR,EK )                                      
!xx      CALL  SUBSPC( B,K,N,IPR,EK,IND,C,C1,C2,OEIC,ESUM,OMEAN,OM )
call  subspc( b,k,n,ek,ind,c,c1,c2,oeic,esum,omean,om )
!                                                                       
!          MODIFICATION OF CROSS-PRODUCTS  X(I,K1) (I=1,K)              
!                                                                       
do 30   i=1,k
!c      B(I) = B(I)*X(I,I)*G(I) / DABS(X(I,I))                            
bb = b(i)*x(i,i)*g(i) / dabs(x(i,i))
d(i) = x(i,k1)
!c   30 X(I,K1) = B(I)                                                    
!xx   30 X(I,K1) = BB
x(i,k1) = bb
30 continue
!                                                                       
!          REGRESSION COEFFICIENTS OF BAYSIAN MODEL                     
!                                                                       
call  recoef( x,k,k,mj1,a )
!                                                                       
do 40   i=1,k
!xx   40 X(I,K1) = D(I)                                                    
x(i,k1) = d(i)
40 continue
!                                                                       
!          RESIDUAL VARIANCE AND AIC                                    
!                                                                       
!c      CALL  SDCOMP( X,A,D,N,K,MJ1,SD )                                  
call  sdcomp( x,a,n,k,mj1,sd )
!                                                                       
if( ipr .eq. 0 )     return
aic = fn*dlog(sd) + 2.d0*ek
!c      WRITE( 6,6 )     SD , EK , AIC                                    
return
!                                                                       
!xx    3 FORMAT( //1H ,18(1H-),/,' BAYESIAN PROCEDURE',/,1H ,18(1H-) )     
!xx    5 FORMAT( 1H ,10D13.5 )                                             
!xx    6 FORMAT( 1H ,'RESIDUAL VARIANCE',16X,'SD =',D19.8,/,1H ,'EQUIVALENT
!xx     1 NUMBER OF PARAMETERS  EK =',F10.3,/,1H ,32X,'AIC =',F15.3 )      
end

