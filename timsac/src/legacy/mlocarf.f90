! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine mlocarf( zs,n,lag,ns0,ksw,nml,zmean,sum,a,mf,sdf,lk0,&
&lk2,sxx,nnf,nns,ms,sdms,aics,mp,sdmp,aicp )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM  MLOCAR                                                   
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
!     TIMSAC 78.3.1.                                                    
!     _                     ___                __                       
!     MINIMUM AIC METHOD OF LOCALLY STATIONARY AR MODEL FITTING; SCALAR 
!                                                                       
!     THIS PROGRAM LOCALLY FITS AUTOREGRESSIVE MODELS TO NON-STATIONARY 
!     SERIES BY MINIMUM AIC PROCEDURE.                                  
!                                                                       
!     BY THIS PROCEDURE, THE DATA OF LENGTH N ARE DIVIDED INTO J LOCALLY
!     STATIONARY SPANS                                                  
!                                                                       
!                <-- N1 --> <-- N2 --> <-- N3 -->          <-- NJ -->   
!               !----------!----------!----------! ...... !----------!  
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
!             REDATA                                                    
!             NONSTA                                                    
!             PRINTA                                                    
!             NRASPE                                                    
!       --------------------------------------------------------------- 
!        INPUTS REQUIRED;                                               
!             MT:       INPUT DEVICE FOR ORIGINAL DATA (MT=5 : CARD READ
!             LAG:      UPPER LIMIT OF THE ORDER OF AR-MODEL, MUST BE LE
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
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4   Z(10000) , TITLE(20)                                   
!c      REAL * 4   TITLE(20)
!c      DIMENSION  X(200,51)                                              
!c      DIMENSION  D(200) , A(50)                                         
!c      DIMENSION  U(51,51)                                               
!C      DIMENSION  TTL(2)                                                 
!C      DATA  TTL / 8H  CURREN,8HT MODEL  /                               
!xx      DIMENSION  ZS(N), Z(N)
!xx      DIMENSION  X(N,LAG+KSW+1)
!xx      DIMENSION  U(LAG+KSW+1,LAG+KSW+1), AA(LAG+KSW), A(LAG+KSW,NML)
!c      REAL*4  TTL(4)                                                 
!c      DATA  TTL / 4H  CU,4HRREN,4HT MO,4HDEL  /                               
!
!xx      DIMENSION  MF(NML), SDF(NML), LK0(NML), LK2(NML)
!xx      DIMENSION  NNF(NML), NNS(NML)
!xx      DIMENSION  MS(NML), SDMS(NML), AICS(NML)
!xx      DIMENSION  MP(NML), SDMP(NML), AICP(NML)
!xx      DIMENSION  SXX(121,NML)
integer n, lag, ns0, ksw, nml, mf(nml), lk0(nml), lk2(nml),&
&nnf(nml), nns(nml), ms(nml), mp(nml)
real(dp) zs(n), zmean, sum, a(lag+ksw,nml), sdf(nml),&
&sxx(121,nml), sdms(nml), aics(nml), sdmp(nml),&
&aicp(nml)
! local
integer i, if, isw, k, l, lk, lk1, mj1, mj2, mx, nf, nj, ns
real(dp) z(n), x(n,lag+ksw+1), u(lag+ksw+1,lag+ksw+1),&
&aa(lag+ksw), b
!                                                                       
!          EXTERNAL SUBROUTINE DECLARATION                              
!                                                                       
external  setx1
!                                                                       
!     PARAMETERS:                                                       
!          ISW:   =0  TO PRODUCE THE MAICE MODEL ONLY (OUTPUTS SUPRRESSE
!                 =1  TO PRODUCE THE MAICE MODEL ONLY                   
!                 =2  TO PRODUCE ALL AR-MODELS (UP TO THE ORDER K)      
!          MJ1:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
!          MJ2:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL                
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
!c      MJ1 = 200                                                         
!c      MJ2 = 51                                                          
mj1 = n
mj2 = lag+1
isw = 0
!
mf(1:nml) = 0
a(1:(lag+ksw),1:nml) = 0.0d0
sxx(1:121,1:nml) = 0.0d0
!                                                                       
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      READ( 5,1 )     LAG , NS , KSW                                    
!
!c      WRITE( 6,3 )                                                      
!c      IF( KSW .EQ. 1 )     WRITE( 6,5 )                                 
!c      IF( KSW .NE. 1 )     WRITE( 6,4 )                                 
!c      WRITE( 6,6 )                                                      
!c      WRITE( 6,2 )    LAG , NS , MT                                     
!
!          ---------------------------------------                      
!          ORIGINAL DATA LOADING AND MEAN DELETION                      
!          ---------------------------------------                      
!
!c      CALL  REDATA( Z,N,MT,TITLE )
call  redata( zs,z,n,zmean,sum )
!c      CLOSE( MT )
!
l  = 0
k  = lag + ksw
mx = k * 2
if = 0
nj = 0
nf = 0
ns = ns0
!                                                                       
100 continue
!                                                                       
lk  = l + k
lk1 = lk + 1
if( lk1     .ge. n  )   go to 200
if( n-lk1   .lt. ns )   ns = n - lk
if( n-lk1-ns.lt. mx )   ns = n - lk
!                                                                       
!          -----------------------------------                          
!          LOCALLY STATIONARY AR-MODEL FITTING                          
!          -----------------------------------                          
!                                                                       
!c      CALL  NONSTA( SETX1,Z,X,U,D,LAG,L,NS,K,IF,ISW,TITLE,MJ1,MJ2,A,MF, 
!c     1SDF )                                                             
nj = nj+1
if ( nj .gt. 1 )  mf(nj) = mf(nj-1)
if ( nj .gt. 1 )  sdf(nj) = sdf(nj-1)
call  nonsta( setx1,z,x,u,lag,l,nf,ns,k,if,isw,mj1,mj2,aa,mf(nj),&
&sdf(nj),nnf(nj),nns(nj),ms(nj),sdms(nj),aics(nj),mp(nj),sdmp(nj),&
&aicp(nj) )
!                                                                       
l = l + ns
!c      IF( IF .EQ. 2 )     LK0 = LK1                                     
if( if .eq. 2 )     lk0(nj) = lk1
if( if .ne. 2 )     lk0(nj) = lk0(nj-1)
!
!          -----------------------                                      
!          PRINT OUT CURRENT MODEL                                      
!          -----------------------                                      
!
!c      LK2 = LK + NS                                                     
!c      CALL  PRINTA( A,SDF,MF,TTL,4,TITLE,LK0,LK2 )                      
lk2(nj) = lk + ns
!                                                                       
!                                                                       
!          ----------------                                             
!          SPECTRUM DISPLAY                                             
!          ----------------                                             
!                                                                       
!c      CALL  NRASPE( SDF,A,B,MF,0,121,TITLE )                            
call  nraspe( sdf(nj),aa,[b],mf(nj),0,120,sxx(:,nj) )
do 110 i = 1,mf(nj)
!xx  110 A(I,NJ) = AA(I)
a(i,nj) = aa(i)
110 continue
!                                                                       
go to 100
!                                                                       
200 continue
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
!
!xx  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//5X,100A)
!xx  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//5X,100A)
!
!c  999 CONTINUE
return
!                                                                       
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( ///1H ,'  FITTING UP TO THE ORDER  K =',I3,'  IS TRIED',/,
!xx     1'   BASIC LOCAL SPAN  NS =',I4,/,'   ORIGINAL DATA INPUT DEVICE  M
!xx     2T =',I3 )                                                         
!xx    3 FORMAT( //  ' PROGRAM TIMSAC 78.3.1',/'   LOCALLY STATIONARY AUTOR
!xx     1EGRESSIVE MODEL FITTING;   SCALAR CASE',//,'   < BASIC AUTOREGRESS
!xx     2IVE MODEL >' )                                                    
!xx    4 FORMAT( 1H ,10X,'Z(I) = A(1)*Z(I-1) + A(2)*Z(I-2) + ... + A(M)*Z(I
!xx     1-M) + E(I)' )                                                     
!xx    5 FORMAT( 1H ,10X,'Z(I) = A(1) + A(2)*Z(I-1) + ... + A(M+1)*Z(I-M) +
!xx     1 E(I)' )                                                          
!xx    6 FORMAT( 1H ,2X,'WHERE',/,11X,'M:     ORDER OF THE MODEL',/,11X,'E(
!xx     1I):  GAUSSIAN WHITE NOISE WITH MEAN 0  AND  VARIANCE SD(M).' )    
!                                                                       
end
!c      SUBROUTINE  NONSTA( SETX,Z,X,U,D,LAG,N0,NS,K,IF,ISW,TITLE,MJ1,MJ2,
!c     1A,MF,SDF )                                                        
subroutine  nonsta( setx,z,x,u,lag,n0,nf,ns,k,if,isw,mj1,mj2,a,mf,&
&sdf,nnf,nns,ms,sdms,aics,mp,sdmp,aicp )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     IN THIS SUBROUTINE THE FOLLOWING TWO MODELS ARE COMPARED AND      
!     THE MODEL WITH LESS AIC IS ACCEPTED AS THE CURRENT MODEL.         
!                                                                       
!       MOVING MODEL:     SUCCESSION OF TWO-MODELS INDEPENDENTLY FITTED 
!                         TO THE DIVIDED DATA                           
!              NF:     DATA LENGTH OF THE PRECEDING STATIONARY BLOCK    
!              NS:     DATA LENGTH OF NEW BLOCK (= BASIC LOCAL SPAN )   
!              AR(MF,SDF):   MAICE AR-MODEL WITH THE ORDER MF AND INNOVA
!                            VARIANCE SDF FITTED TO THE PRECEDING STATIO
!                            BLOCK                                      
!              AR(MS,SDS):   MAICE AR-MODEL (ORDER MS AND INNOVATION VAR
!                            SDS) FITTED TO THE NEWLY OBTAINED DATA     
!              AICS = NF*LOG(SDF) + NS*LOG(SDS) + 2*(MF+MS+2)           
!                                                                       
!                                AR(MF,SDF)       AR(MS,SDS)            
!                             !---------------!----------------!        
!                              <---- NF -----> <----- NS ----->         
!                                                                       
!                                                                       
!       CONSTANT MODEL:   AR MODEL FITTED TO THE POOLED DATA            
!              NP:     DATA LENGTH OF POOLED DATA  (=NF+NS)             
!              AR(MP,SDP):   MAICE AR-MODEL FITTED TO THE POOLED DATA   
!                            (ORDER = MP, INNOVATION VARIANCE = SDP)    
!              AICP = NP*LOG(SDP) + 2*(MP+1)                            
!                                                                       
!                                        AR(MP,SDP)                     
!                             !--------------------------------!        
!                              <------------ NP -------------->         
!                                                                       
!                                                                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             COPY                                                      
!             ARMFIT                                                    
!             HUSHLD                                                    
!             REDUCT                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          SETX:   EXTERNAL SUBROUTINE DESIGNATION                      
!          Z:      ORIGINAL DATA VECTOR                                 
!          X:      WORKING AREA                                         
!          U:      WORKING AREA                                         
!          D:      WORKING AREA                                         
!          LAG:    UPPER LIMIT OF THE ORDER OF AR-MODEL                 
!          NS:     LENGTH OF BASIC LOCAL SPAN                           
!          IF:     =0   FIT INITIAL AR-MODEL AND STORE                  
!                  >0   UPDATE THE CURRENT MODEL                        
!          ISW:    =0  TO PRODUCE THE MAICE MODEL ONLY (OUTPUTS SUPPRESS
!                  =1  TO PRODUCE THE MAICE MODEL ONLY                  
!                  =2  TO PRODUCE ALL AR-MODELS (UP TO THE ORDER K)     
!          TITLE:  TITLE OF DATA                                        
!          MJ1:    ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM          
!          MJ2:    ABSOLUTE DIMENSION OF U IN THE MAIN PROGRAM          
!                                                                       
!       OUTPUTS:                                                        
!          A:      AR-COEFFICIENTS OF THE CURRENT MODEL                 
!          MF:     ORDER OF THE CURRENT MODEL                           
!          SDF:    INNOVATION VARIANCE OF THE CURRENT MODEL             
!          IF:     =1   MODEL UNSWITCHED                                
!                  =2   MODEL SWITCHED                                  
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4  Z(1) , TITLE(1)                                         
!c      REAL * 4   TITLE(1)
!c      DIMENSION  X(MJ1,1) , U(MJ2,1) , A(1)                      
!c      DIMENSION  B(50)                                                  
!x      DIMENSION  Z(1)
!xx      DIMENSION  Z(MJ1)
!xx      DIMENSION  X(MJ1,K+1) , U(K+1,K+1) , A(K)                      
!xx      DIMENSION  B(K)                                                  
!xx      DIMENSION  SDS(K+1), AS(K+1), DICS(K+1)
!xx      DIMENSION  SDP(K+1), AP(K+1), DICP(K+1)
integer lag, n0, nf, ns, k, if, isw, mj1, mj2, mf, nnf, nns,&
&ms, mp
real(dp) z(mj1), x(mj1,k+1), u(k+1,k+1), a(k), sdf, sdms,&
&aics, sdmp, aicp
! local
integer i, k1, k2, np
real(dp) b(k), sds(k+1), as(k+1), dics(k+1), sdp(k+1),&
&ap(k+1), dicp(k+1), aicms, aicmp
external  setx
!                                                                       
k1 = k + 1
k2 = k1*2
nnf = 0
nns = 0
!x      IFG = 0
!                                                                     +-
!       ---  DATA LOADING AND HOUSEHOLDER TRANSFORMATION  ---         ! 
!                                                                     +-
!c      CALL  REDUCT( SETX,Z,D,NS,N0,K,MJ1,LAG,X )                        
call  reduct( setx,z,ns,n0,k,mj1,lag,x )
!                                                                     +-
!       ---  AR-MODEL FITTING TO NEW SET OF DATA  ---                 ! 
!                                                                     +-
!c      CALL  ARMFIT( X,K,LAG,NS,ISW,TITLE,MJ1,B,SDS,MS )                 
!x      CALL  ARMFIT( X,K,LAG,NS,ISW,MJ1,B,MS,SDS,AS,DICS,SDMS,AICMS,
!x     *              IFG,LU )
call  armfit( x,k,lag,ns,isw,mj1,b,ms,sds,as,dics,sdms,aicms)
!                                                                NO  +--
!                                                              +-----!  
!                                                              !     +--
if( if .ne. 0 )     go to 10
!                                                              !        
!                                                              !      +-
!       ---  MAKE A COPY OF X ON U  ---                        !      ! 
!                                                              !      +-
call  copy( x,k1,0,0,mj1,mj2,u )
!                                                              !        
!       ---  AIC FOR INITIAL LOCAL MODEL  ---                  !        
!                                                              !        
!c      AICS = NS*DLOG(SDS) + 2.D0*(MS+1)                                 
aics = ns*dlog(sdms) + 2.d0*(ms+1)
!                                                              !    ( GO
!c      WRITE( 6,5 )     NS , MS , SDS , AICS                             
nns = ns
go to 20
!                                                              !        
!                                                              +--------
!       ---  AIC FOR MOVING MODEL  ---                                  
!                                                                       
!c   10 AICS = NF*DLOG(SDF) + NS*DLOG(SDS) + 2.D0*(MF+MS+2)               
10 aics = nf*dlog(sdf) + ns*dlog(sdms) + 2.d0*(mf+ms+2)
!                                                                       
!c      WRITE( 6,4 )                                                      
!c      WRITE( 6,6 )     NF , NS , MS , SDS , AICS                        
nnf = nf
nns = ns
!                                                                       
!          -------------------------------                              
!          AR-MODEL FITTING TO POOLED DATA                              
!          -------------------------------                            +-
!                                                                     ! 
!       ---  MAKE COPIES OF U AND X  ---                              ! 
!                                                                     +-
call  copy( x,k1,0,k2,mj1,mj1,x )
call  copy( u,k1,0,k1,mj2,mj1,x )
!                                                                     +-
!       HOUSEHOLDER TRANSFORMATION  ---                               ! 
!                                                                     +-
!c      CALL  HUSHLD( X,D,MJ1,K2,K1 )                                     
call  hushld( x,mj1,k2,k1 )
!                                                                       
!       ---  AR MODEL FITTING TO POOLED DATA  ---                       
!                                                                       
np = nf + ns
!c      CALL  ARMFIT( X,K,LAG,NP,ISW,TITLE,MJ1,A,SDP,MP )                 
!x      CALL  ARMFIT( X,K,LAG,NP,ISW,MJ1,A,MP,SDP,AP,DICP,SDMP,AICMP,
!x     *              IFG,LU )
call  armfit( x,k,lag,np,isw,mj1,a,mp,sdp,ap,dicp,sdmp,aicmp)
!                                                                       
!       ---  AIC FOR CONSTANT MODEL  ---                                
!                                                                       
!c     AICP = NP*DLOG(SDP) + 2.D0*(MP+1)                                 
!c      WRITE( 6,7 )     NP , MP , SDP , AICP                             
aicp = np*dlog(sdmp) + 2.d0*(mp+1)
!                                                                       
!          --------------------                              YES  +-----
!          COMPARISON OF MODELS                             +-----!AICS 
!          --------------------                             !     +-----
!                                                           !           
if( aics .ge. aicp )     go to  40
!                                                           !           
!          ------------------------                                     
!          SPECTRUM CHANGE DETECTED                                     
!          ------------------------                                     
!                                                           !           
!c      WRITE( 6,8 )                                                      
!                                                           !         +-
!       ---  MAKE A COPY OF X2 ON U  ---                    !         !X
!                                                           !         +-
call  copy( x,k1,k2,0,mj1,mj2,u )
!                                                           !           
20 if = 2
nf = ns
mf = ms
!                                                           !           
do 30  i=1,mf
!xx   30 A(I) = B(I)
a(i) = b(i)
30 continue
!                                                           !           
!c      SDF = SDS                                                         
sdf = sdms
go to 50
!                                                           +-----------
40 if = 1
!                                                                       
!          ------------------                                   (DATA PO
!          SPECTRUM UNCHANGED                                           
!          ------------------                                           
!                                                                   +---
!       ---  MAKE A COPY OF X ON U  ---                             ! X-
!                                                                   +---
!                                                                       
call  copy( x,k1,0,0,mj1,mj2,u )
!                                                                       
!c      WRITE( 6,9 )                                                      
!c      SDF = SDP                                                         
sdf = sdmp
mf = mp
nf = nf + ns
!                                                                       
50 continue
!
return
!
!xx    4 FORMAT( //1H ,'---  THE FOLLOWING TWO MODELS ARE COMPARED  ---' ) 
!xx    5 FORMAT( //1H ,'INITIAL LOCAL MODEL:    NS =',I5,5X,'MS =',I3,5X,  
!xx     1  'SDS =',D16.8,5X,'AICS =',F16.3 )                               
!xx    6 FORMAT( 1H ,'MOVING MODEL:      (NF =',I5,', NS =',I4,1H),5X,'MS =
!xx     2',I3,5X,'SDS =',D16.8,5X,'AICS =',F16.3 )                         
!xx    7 FORMAT( 1H ,'CONSTANT MODEL:    (NP =',I5,1H),15X,'MP =',I3,5X,'SD
!xx     3P =',D16.8,5X,'AICP =',F16.3 )                                    
!xx    8 FORMAT( //1H ,37(1H*),/,1H ,'*****',27X,'*****',/,1H ,'*****     N
!xx     1EW MODEL ADOPTED     *****',/,1H ,'*****',27X,'*****',/,1H ,37(1H*
!xx     2) )                                                               
!xx    9 FORMAT( 1H ,'*****  CONSTANT MODEL ADOPTED  *****' )              
!xx  700 FORMAT( 1H ,'-----  X  -----' )                                   
!xx  620 FORMAT( 1H ,10D13.5 )                                             
!xx  600 FORMAT( 1H ,'N =',I5,5X,'K =',I5,5X,'M =',I5,5X,'MT =',I5,5X,     
!xx     * 'DATA FORMAT =',15A4 )                                           
!xx  601 FORMAT( 1H ,'-----  ORIGINAL DATA  -----',/,(1X,10D13.5) )        
!                                                                       
end
