! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine  mlomarf( zs,n,id,c,lag,ns0,ksw,k,zmean,zvari,nf,ns,ms,&
&aic,mp,aicp,mf,aicf,a,e,lk0,lke,m )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM  MLOMAR                                                   
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
!     TIMSAC 78.3.3.                                                    
!     _                     __                 _            __          
!     MINIMUM AIC METHOD OF LOCALLY STATIONARY MULTIVARIATE AR MODEL FIT
!                                                                       
!     THIS PROGRAM LOCALLY FITS MULTI-VARIATE AUTOREGRESSIVE MODELS TO  
!     NON-STATIONARY TIME SERIES BY THE MINIMUM AIC PROCEDURE USING THE 
!     HOUSEHOLDER TRANSFORMATION.                                       
!                                                                       
!     BY THIS PROCEDURE, THE DATA OF LENGTH N ARE DIVIDED INTO J LOCALLY
!     STATIONARY SPANS                                                  
!                                                                       
!                <-- N1 --> <-- N2 --> <-- N3 -->          <-- NJ -->   
!               !----------!----------!----------!--------!----------!  
!                <-----------------------  N  ---------------------->   
!                                                                       
!     WHERE NI (I=1,...,J) DENOTES THE NUMBER OF BASIC SPANS, EACH OF   
!     LENGTH NS, WHICH CONSTITUTE THE I-TH LOCALLY STATIONARY SPAN.     
!     AT EACH LOCAL SPAN, THE PROCESS IS REPRESENTED BY A STATIONARY    
!     AUTOREGRESSIVE MODEL.                                             
!                                                                       
!                                                                       
!       --------------------------------------------------------------- 
!       REFERENCE:                                                      
!          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
!          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
!          30,B,351-363.                                                
!       --------------------------------------------------------------- 
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
!             MRDATA                                                    
!             MNONST                                                    
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
!c      DIMENSION  X(200,100) , U(100,100) , D(200)                       
!c      DIMENSION  A(5,5,50) , B(5,5,50) , E(5,5)                         
!xx      DIMENSION  ZS(N,ID), Z(N,ID), C(ID)
!xx      DIMENSION  ZMEAN(ID), ZVARI(ID)
!xx      DIMENSION  A(ID,ID,LAG,K) , B(ID,ID,LAG) , E(ID,ID,K)             
!xx      DIMENSION  NF(K), NS(K), MS(K), MP(K), MF(K)
!xx      DIMENSION  AIC(K), AICP(K), AICF(K)
!xx      DIMENSION  LK0(K), LKE(K)
!xxcx      DIMENSION  X(N,((LAG+1)*ID+KSW)*2)
!xx      DIMENSION  X(((LAG+1)*ID+KSW)*4,((LAG+1)*ID+KSW)*2)
!xx      DIMENSION  U(((LAG+1)*ID+KSW)*2,((LAG+1)*ID+KSW)*2)
integer n, id, lag, ns0, ksw, k, nf(k), ns(k), ms(k), mp(k),&
&mf(k), lk0(k), lke(k), m
real(dp) zs(n,id), c(id), zmean(id), zvari(id), aic(k),&
&aicp(k), aicf(k), a(id,id,lag,k), e(id,id,k)
! local
integer if, kd, l, lk, lk1, mj1, mj2, mj3, mx, nnf
real(dp) z(n,id), b(id,id,lag),&
&x(((lag+1)*id+ksw)*4,((lag+1)*id+ksw)*2),&
&u(((lag+1)*id+ksw)*2,((lag+1)*id+ksw)*2)
!
!c      CHARACTER(100)  IFLNAM,OFLNAM
!c      CALL FLNAM2( IFLNAM,OFLNAM,NFL )
!c      IF ( NFL.EQ.0 ) GO TO 999
!c      IF ( NFL.EQ.2 ) THEN
!c         OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR )
!c      ELSE
!c         CALL SETWND
!c      END IF
!                                                                       
!       PARAMETERS:                                                     
!          MJ:    ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
!          MJ1:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
!          MJ2:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
!          MJ3:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
!                                                                       
!c      MJ = 1500                                                         
!c      MJ1 = 200                                                         
!c      MJ2 = 100                                                         
!c      MJ3 = 5                                                           
mj2 = ((lag+1)*id+ksw)*2
mj1 = mj2*2
mj3 = id
!
nf(1:k) = 0
ns(1:k) = 0
ms(1:k) = 0
aic(1:k) = 0.0d0
mp(1:k) = 0
aicp(1:k) = 0.0d0
mf(1:k) = 0
aicf(1:k) = 0.0d0
a(1:id,1:id,1:lag,1:k) = 0.0d0
e(1:id,1:id,1:k) = 0.0d0
lk0(1:k) = 0
lke(1:k) = 0
x(1:mj1,1:mj2) = 0.0d0
u(1:mj2,1:mj2) = 0.0d0
!
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      READ( 5,1 )     LAG , NS , KSW                                    
!                                                                       
ns(1) = ns0
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
!                                                                       
!                                                                       
if = 0
m = 0
nf(1) = 0
111 continue
!xx      M = M+1
lk = l + lag
lk1 = lk + 1
if( lk1 .ge. n )     go to 300
m = m+1
!c      IF( N-LK1 .LE. NS )     NS = N - LK                               
!c      IF( N-LK1-NS .LT. MX )     NS = N - LK                            
if( m.ne. 1 )  then
aicf(m) = aicf(m-1)
ns(m) = ns(m-1)
lk0(m) = lk0(m-1)
end if
if( n-lk1 .le. ns(m) )     ns(m) = n - lk
if( n-lk1-ns(m) .lt. mx )     ns(m) = n - lk
!                                                                       
!c      CALL  MNONST( Z,X,U,D,KSW,LAG,L,NS,ID,IF,MJ,MJ1,MJ2,MJ3,A,B,E,MF, 
!c     *              AIC )                                               
call mnonst( z,x,u,ksw,lag,l,nnf,nf(m),ns(m),id,if,n,mj1,mj2,mj3,&
&a(1,1,1,m),b,e(1,1,m),ms(m),aic(m),mp(m),aicp(m),mf(m),aicf(m) )
!                                                                      
!c      L = L + NS                                                        
!c      IF( IF .EQ. 2 )     LK0 = LK1                                     
l = l + ns(m)
if( if .eq. 2 )     lk0(m) = lk1
!                                                                       
!c      LKE = LK + NS                                                     
lke(m) = lk + ns(m)
!c      WRITE( 6,13 )                                                     
!c      WRITE( 6,16 )                                                     
!c      WRITE( 6,14 )     LK0 , LKE                                       
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
!c  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
!c      GO TO 999
!
!c  910 CONTINUE
!c      IF ( NFL.EQ.2 ) CLOSE( 6 )
!c#ifdef __linux__
!cC	reopen #6 as stdout
!c      IF ( NFL.EQ.2 ) OPEN(6, FILE='/dev/fd/1')
!c#endif
!cC /* __linux__ */
!c      WRITE(6,610) IVAR,IFLNAM
!c  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
!
!c  999 CONTINUE
!c      close(3)
return
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( ///1H ,'PROGRAM TIMSAC 78.3.3',/,'   LOCALLY STATIONARY MU
!xx     1LTI-VARIATE AUTOREGRESSIVE MODEL FITTING;',//,'  < BASIC AUTOREGRE
!xx     2SSIVE MODEL >' )                                                  
!xx    3 FORMAT( ///1H ,'  FITTING UP TO THE ORDER  K =',I3,'  IS TRIED',/,
!xx     1'   BASIC LOCAL SPAN LENGTH  NS =',I4,/,'   ORIGINAL DATA INPUT DE
!xx     2VICE  MT =',I3 )                                                  
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
!c      SUBROUTINE  MNONST( Z,X,U,D,KSW,LAG,N0,NS,ID,IF,MJ,MJ1,MJ2,MJ3,   
!c     1                    A,B,E,MF,AICF )                               
subroutine  mnonst( z,x,u,ksw,lag,n0,nnf,nf,ns,id,if,mj,mj1,mj2,&
&mj3,a,b,e,ms,aicfs,mp,aicp,mf,aicf )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     IN THIS SUBROUTINE THE FOLLOWING TWO MODELS ARE COMPARED AND      
!     THE MODEL WITH THE SMALLER AIC IS ACCEPTED AS THE CURRENT MODEL.  
!                                                                       
!       MOVING MODEL:      SUCCESSION OF TWO AR-MODELS INDEPENDENTLY FIT
!                          TO THE FORMER AND PRESENT BLOCK OF DATA      
!          NF:    DATA LENGTH OF THE PRECEDING STATIONARY BLOCK         
!          NS:    DATA LENGTH OF NEW BLOCK (A BASIC LOCAL SPAN)         
!          AR(MF,SDF):  MAICE AR-MODEL WITH THE ORDER MF AND INNOVATION 
!                       VARIANCE SDF FITTED TO THE PRECEDING STATIONARY 
!                       BLOCK                                           
!          AR(MS,SDS):  MAICE AR-MODEL (ORDER MS AND INNOVATION VARIANCE
!                       SDS) FITTED TO THE NEWLY OBTAINED DATA          
!                                                                       
!                                AR(MF,SDF)      AR(MS,SDS)             
!                             !---------------!---------------!         
!                              <---- NF -----> <----- NS ---->          
!                                                                       
!       CONSTANT MODEL:    AR MODEL FITTED TO THE POOLED DATA OF THE FOR
!                          AN PRESENT BLOCK                             
!          NP:    DATA LENGTH OF POOLED DATA                            
!          AR(MP,SDP):  MAICE AR-MODEL FITTED TO THE POOLED DATA        
!                       (ORDER=MP, INNOVATION VARIANCE=SDP)             
!                                                                       
!                                       AR(MP,SDP)                      
!                             !-------------------------------!         
!                              <----------- NP -------------->          
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             COPY                                                      
!             HUSHLD                                                    
!             MARFIT                                                    
!             MREDCT                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          Z:      ORIGINAL DATA; Z(K,I) (K=1,N) REPRESENTS THE RECORD O
!                  THE I-TH CHANNEL                                     
!          X:      WORKING AREA                                         
!          U:      WORKING AREA                                         
!          D:      WORKING AREA                                         
!          KSW:    =0   CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR  
!                  =1   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESS
!          LAG:    UPPER LIMIT OF THE ORDER OF AR MODEL                 
!          N0:     INDEX OF THE END POINT OF THE FORMER SPAN            
!          NS:     LENGTH OF BASIC LOCAL SPAN                           
!          ID:     DIMENSION OF DATA                                    
!          IF:                                                          
!          MJ:     ABSOLUTE DIMENSION OF Z IN THE MAIN PROGRAM          
!          MJ1:    ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM          
!          MJ2:    ABSOLUTE DIMENSION OF U IN THE MAIN PROGRAM          
!          MJ3:    ABSOLUTE DIMENSION OF A IN THE MAIN PROGRAM          
!          B:      WORKING AREA                                         
!                                                                       
!       OUTPUTS:                                                        
!          A:      AR COEFFICIENT MATRICES OF THE CURRENT MODEL         
!          MF:     ORDER OF THE CURRENT MODEL                           
!          AICF:   AIC OF THE CURRENT MODEL                             
!                                                                       
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4  Z                                                       
!c      DIMENSION  X(MJ1,1) , U(MJ2,1) , D(1)                             
!x      DIMENSION  Z(MJ,1)                                                
!xx      DIMENSION  Z(MJ,ID)                                                
!xx      DIMENSION  X(MJ1,MJ2) , U(MJ2,MJ2)
!xx      DIMENSION  A(ID,ID,LAG) , B(ID,ID,LAG) , E(ID,ID)                 
!xx      DIMENSION  AI(ID,ID,LAG), BI(ID,ID,LAG), EI(ID,ID) 
!c      DIMENSION  Y(100,100)                                             
!c      DIMENSION  C(10) , EX(10)                                         
!xx      DIMENSION  C(ID) , EX(ID)
!xx      DIMENSION  AIC(LAG+1,ID), SD(LAG+1,ID), DIC(LAG+1,ID)
!xx      DIMENSION  AICM(ID), SDM(ID), M(ID)
!xx      DIMENSION  JNDF(MJ2,ID), AF(MJ2,ID), NPR(ID), AAIC(ID)
integer ksw, lag, n0, nnf, nf, ns, id, if, mj, mj1, mj2, mj3,&
&ms, mp, mf
real(dp) z(mj,id), x(mj1,mj2), u(mj2,mj2), a(id,id,lag),&
&b(id,id,lag), e(id,id), aicfs, aicp, aicf
! local
integer i, ii, ipr, j, k1, kd1, kd2, mj4, np, m(id), jndf(mj2,id),&
&npr(id)
real(dp) ai(id,id,lag), bi(id,id,lag), ei(id,id), c(id),&
&ex(id), aic(lag+1,id), sd(lag+1,id),&
&dic(lag+1,id), aicm(id), sdm(id), af(mj2,id),&
&aaic(id), aics
!                                                                       
!                                                                       
!c      MJ4 = 50                                                          
mj4 = lag
ipr = 0
!x      IFG = 0
!                                                                       
k1 = lag + 1
kd1 = k1*id + ksw
kd2 = kd1 * 2
!                                                                       
!       HOUSEHOLDER'S REDUCTION                                         
!                                                                       
!c      CALL  MREDCT( Z,D,NS,N0,LAG,ID,MJ,MJ1,KSW,X )                     
call  mredct( z,ns,n0,lag,id,mj,mj1,ksw,x )
!                                                                       
!                                                                       
!       AR-MODEL FITTING BY THE MINIMUM AIC PROCEDURE                   
!                                                                       
!c      CALL  MARFIT( X,Y,D,NS,ID,LAG,KSW,MJ1,MJ3,MJ4,MJ2,0,IPR,B,E,EX,C, 
!c     *             MS,AICS )                                            
!xx      CALL  MARFIT( X,NS,ID,LAG,KSW,MJ1,MJ3,MJ4,MJ2,0,IPR,AIC,SD,DIC,
call  marfit( x,ns,id,lag,ksw,mj1,mj3,mj4,kd1,0,ipr,aic,sd,dic,&
!x     *AICM,SDM,M,BI,EI,B,E,EX,C,MS,AICS,JNDF,AF,NPR,AAIC,IFG,LU )
&aicm,sdm,m,bi,ei,b,e,ex,c,ms,aics,jndf,af,npr,aaic )
!                                                                       
if( if .ne. 0 )     go to 10
!                                                                       
call  copy( x,kd1,0,0,mj1,mj2,u )
!                                                                       
!c      WRITE( 6,5 )     NS , MS , AICS                                   
go to 20
!                                                                       
!                                                                       
!c   10 AIC = AICF + AICS                                                 
10 aicfs = aicf + aics
!                                                                       
!c      WRITE( 6,4 )                                                      
!c      WRITE( 6,6 )     NF , NS , MS , AIC                               
nf = nnf
!                                                                       
call  copy( x,kd1,0,kd2,mj1,mj1,x )
call  copy( u,kd1,0,kd1,mj2,mj1,x )
!                                                                       
!       ---  HOUSEHOLDER TRANSFORMATION  ---                            
!                                                                       
!c      CALL  HUSHLD( X,D,MJ1,KD2,KD1 )                                   
call  hushld( x,mj1,kd2,kd1 )
!                                                                       
!                                                                       
!                                                                       
!       ---  AR-MODEL FITTING FOR POOLED DATA  ---                      
!                                                                       
np = nnf + ns
!c      CALL  MARFIT( X,Y,D,NP,ID,LAG,KSW,MJ1,MJ3,MJ4,MJ2,0,IPR,A,E,EX,C, 
!c     *              MP,AICP )                                           
!xx      CALL  MARFIT( X,NP,ID,LAG,KSW,MJ1,MJ3,MJ4,MJ2,0,IPR,AIC,SD,DIC, 
call  marfit( x,np,id,lag,ksw,mj1,mj3,mj4,kd1,0,ipr,aic,sd,dic,&
!x     *AICM,SDM,M,AI,EI,A,E,EX,C,MP,AICP,JNDF,AF,NPR,AAIC,IFG,LU )
&aicm,sdm,m,ai,ei,a,e,ex,c,mp,aicp,jndf,af,npr,aaic )
!                                                                       
!c      WRITE( 6,7 )     NP , MP , AICP                                   
!                                                                       
!c      IF( AIC .GE. AICP )     GO TO 40                                  
if( aicfs .ge. aicp )     go to 40
!                                                                       
!c      WRITE( 6,8 )                                                      
!                                                                       
call  copy( x,kd1,kd2,0,mj1,mj2,u )
!                                                                       
20 continue
if = 2
nnf = ns
mf = ms
aicf = aics
!                                                                       
!xx      DO 30  II=1,MF
!xx      DO 30  J=1,ID
do 32  ii=1,mf
do 31  j=1,id
do 30  i=1,id
!xx   30 A(I,J,II) = B(I,J,II) 
a(i,j,ii) = b(i,j,ii)
30 continue
31 continue
32 continue
go to 50
!                                                                       
!                                                                       
!                                                                       
40 continue
if = 1
call  copy( x,kd1,0,0,mj1,mj2,u )
!                                                                       
!c      WRITE( 6,9 )                                                      
nnf = nnf + ns
mf = mp
aicf = aicp
!                                                                       
!                                                                       
50 continue
!                                                                       
!                                                                       
return
!                                                                       
!xx  600 FORMAT( 1H ,'N =',I5,5X,'ID =',I5,5X,'K =',I5,5X,'M =',I5,5X,     
!xx     U  'MT =',I5,5X,'DATA FORMAT =',10A4 )                             
!xx  601 FORMAT( 1H ,'-----  ORIGINAL DATA  -----' )                       
!xx  602 FORMAT( 1H ,9X,'I',10X,'MEAN',7X,'VARIANCE' )                     
!xx  610 FORMAT( 1H ,10D13.5 )                                             
!xx  620 FORMAT( 1H ,I10,2D15.7 )                                          
!xx    2 FORMAT( 20A4 )                                                    
!xx    3 FORMAT( 8F10.0 )                                                  
!xx    4 FORMAT( //1H ,'---  THE FOLLOWING TWO MODELS ARE COMPARED  ---' ) 
!xx    5 FORMAT( //1H ,'INITIAL LOCAL MODEL:   NS =',I5,5X,'MS =',I3,5X,   
!xx     1'AIC =',F16.3 )                                                   
!xx    6 FORMAT( 1H ,'MOVING MODEL:     (NF =',I5,', NS =',I4,')',5X,      
!xx     1 'MS =',I3,5X,'AIC =',F16.3 )                                     
!xx    7 FORMAT( 1H ,'CONSTANT MODEL:   (NP =',I5,')',15X,'MP =',I3,5X,'AIC
!xx     1 =',F16.3 )                                                       
!xx    8 FORMAT( //1H ,37(1H*),/,1H ,'*****',27X,'*****',/,1H ,'*****     N
!xx     1EW MODEL ADOPTED     *****',/,1H ,'*****',27X,'*****',/,1H ,37(1H*
!xx     2) )                                                               
!xx    9 FORMAT( 1H ,'*****  CONSTANT MODEL ADOPTED  *****' )              
!xx   19 FORMAT( 1H ,'*****',27X,'*****' )                                 
!xx   21 FORMAT( 1H ,37(1H*) )                                             
!xx   22 FORMAT( 1H ,// )                                                  
!xx   24 FORMAT( 1H ,'LK1 =',I5,5X,'M =',I5,/,1H ,130(1H*) )               
!                                                                       
end
