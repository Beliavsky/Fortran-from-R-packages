! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine unibarf( zs,n,lag,zmean,sum,sd,aic,dic,imin,aicm,sdmin,&
&b1,c,d,b2,aicb,sdb,pn,a,sxx )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM  UNIBAR                                                   
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
!     TIMSAC 78.1.2.                                                    
!     ___        _                  __                                  
!     UNIVARIATE BAYESIAN METHOD OF AR MODEL FITTING                    
!                                                                       
!     THIS PROGRAM FITS AN AUTOREGRESSIVE MODEL BY A BAYESIAN PROCEDURE.
!     THE LEAST SQUARES ESTIMATES OF THE PARAMETERS ARE OBTAINED BY THE 
!     HOUSEHOLDER TRANSFORMATION.                                       
!     ----------------------------------------------------------------- 
!     THE BASIC STATISTIC AIC IS DEFINED BY                             
!                                                                       
!               AIC = N * LOG( SD )  +  2 * M  ,                        
!     WHERE                                                             
!       N:     DATA LENGTH,                                             
!       SD:    ESTIMATE OF INNOVATION VARIANCE, THE AVERAGE OF THE      
!              SQUARED RESIDUALS                                        
!       M:     ORDER OF THE MODEL                                       
!     ----------------------------------------------------------------- 
!     BAYESIAN WEIGHT OF THE M-TH ORDER MODEL IS DEFINED BY             
!            W(M)  =  CONST  *  C(M)  /  ( M + 1 ),                     
!     WHERE                                                             
!            CONST =  NORMALIZING CONSTANT                              
!            C(M)  =  EXP( -0.5*AIC(M) ).                               
!     ----------------------------------------------------------------- 
!     THE EQUIVALENT NUMBER OF FREE PARAMETERS FOR THE BAYESIAN MODEL IS
!     DEFINED BY                                                        
!            EK    =  D(1)**2 + ... + D(K)**2 + 1                       
!     WHERE D(J) IS DEFINED BY                                          
!            D(J) = W(J) + ... + W(K).                                  
!                                                                       
!     M IN THE DEFINITION OF AIC IS REPLACED BY EK TO DEFINE AN EQUIVALE
!     FOR A BAYESIAN MODEL.                                             
!     ----------------------------------------------------------------- 
!       REFERENCES:                                                     
!          H.AKAIKE(1978), "A BAYESIAN EXTENSION OF THE MINIMUM AIC     
!          PROCEDURE OF AUTOREGRESSIVE MODEL FITTING.",  RESEARCH MEMO. 
!          NO. 126, THE INSTITUTE OF STATISTICAL MATHEMATICS; TOKYO.   
!          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
!          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
!          30,B,351-363.                                                
!       --------------------------------------------------------------- 
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
!             REDATA                                                    
!             REDUCT                                                    
!             ARBAYS                                                    
!             NRASPE                                                    
!       --------------------------------------------------------------- 
!       INPUTS REQUIRED:                                                
!             LAG         :  ORDER OF THE AR-MODEL, MUST BE LESS THAN 10
!             MT          :  INPUT DEVICE FOR ORIGINAL DATA (MT=5 : CARD
!                                                                       
!               -- THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDA
!                                                                       
!             TITLE:    SPECIFICATION OF DATA                           
!             N:        DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000
!             DFORM:    INPUT DATA FORMAT SPECIFICATION STATEMENT.      
!                       -- FOR EXAMPLE --     (8F10.5)                  
!             (Z(I),I=1,N):  ORIGINAL DATA                              
!               --------------------------------------------------------
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4   Z(10000) , TITLE(20)                                   
!c      REAL * 4   TITLE(20)
!c      DIMENSION  Z(10000)
!c      DIMENSION  X(200,101) , D(200) , A(100) , B(100)                  
!xx      DIMENSION  ZS(N), Z(N)
!xx      DIMENSION  X(N-LAG,LAG+1), D(LAG), A(LAG), B1(LAG), B2(LAG)
!xx      DIMENSION  SD(LAG+1), AIC(LAG+1), DIC(LAG+1), C(LAG+1)
!xx      DIMENSION  SXX(121)
integer n, lag, imin
real(dp) zs(n), zmean, sum, sd(lag+1), aic(lag+1),&
&dic(lag+1), aicm, sdmin, b1(lag), c(lag+1),&
&d(lag), b2(lag), aicb, sdb, pn, a(lag), sxx(121)
! local
integer isw, k, mj1, n1, nmk
real(dp) z(n), x(n-lag,lag+1), b
!                                                                       
!        EXTERNAL SUBROUTINE DECLARATION:                               
!                                                                       
external  setx1
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
!        PARAMETERS:                                                    
!             MJ1:  ABSOLUTE DIMENSION FOR SUBROUTINE CALL              
!             ISW:  =0   OUTPUTS ARE SUPPRESSED                         
!                   >0   OUTPUTS ARE PRINTED OUT                        
!c      MJ1 = 200                                                         
mj1 = n-lag
isw = 1
!                                                                       
!c      WRITE( 6,3 )                                                      
!c      WRITE( 6,4 )                                                      
!
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )                                                                     
!C      READ( 5,1 )     MT                                                
!c      READ( 5,1 )     LAG                                               
!c      WRITE( 6,2 )    LAG , MT                                          
!                                                                       
!                                                                       
!          +---------------------------------------+                  +-
!          ! ORIGINAL DATA INPUT AND MEAN DELETION !                  ! 
!          +---------------------------------------+                  +-
!
!c      CALL  REDATA( Z,N,MT,TITLE )                                      
!c      CLOSE( MT )
call redata( zs,z,n,zmean,sum )
k = lag
nmk = n-lag
!                                                                       
!          +-----------------------+                                  +-
!          ! HOUSEHOLDER REDUCTION !                                  ! 
!          +-----------------------+                                  +-
!                                                                       
n1 = n+1
!c      CALL  REDUCT( SETX1,Z,D,NMK,0,K,MJ1,LAG,X )                       
call  reduct( setx1,z,nmk,0,k,mj1,lag,x )
!                                                                       
!         +-----------------------------------------------+           +-
!         ! AUTOREGRESSIVE MODEL FITTING (BAYESIAN MODEL) !           ! 
!         +-----------------------------------------------+           +-
!                                                                       
!c      CALL  ARBAYS( X,D,K,LAG,NMK,ISW,TITLE,MJ1,A,B,SDB,AICB )          
!xx      CALL  ARBAYS( X,D,K,LAG,NMK,ISW,MJ1,SD,AIC,DIC,AICM,SDMIN,IMIN,
call  arbays( x,d,k,nmk,isw,mj1,sd,aic,dic,aicm,sdmin,imin,&
&a,b1,b2,c,sdb,pn,aicb )
!                                                                       
!          +------------------------+                                 +-
!          ! POWER SPECTRUM DISPLAY !                                 ! 
!          +------------------------+                                 +-
!                                                                       
!c      CALL  NRASPE( SDB,A,B,K,0,120,TITLE )                             
call  nraspe( sdb,a,[b],k,0,120,sxx )
!c      GO TO 999
!
!c  900 CONTINUE
!c      WRITE(6,600) IVAR,OFLNAM
!c  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
!c       GO TO 999
!
!c  910 CONTINUE
!c      IF ( NFL.EQ.2 ) CLOSE( 6 )
!c#ifdef __linux__
!      reopen #6 as stdout
!c      IF ( NFL.EQ.2 ) OPEN(6, FILE='/dev/fd/1')
!c#endif
! /* __linux__ */
!c      WRITE(6,610) IVAR,IFLNAM
!c  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
!                                                                       
!c  999 CONTINUE
return
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( 1H ,I4,'-TH ORDER BAYESIAN MODEL IS FITTED',/,1H ,2X,     
!xx     1'ORIGINAL DATA INPUT DEVICE  MT =',I4 )                           
!xx    3 FORMAT( ' PROGRAM TIMSAC 78.1.2',/'   EXPONENTIALLY WEIGHTED BAYES
!xx     1IAN AUTOREGRESSIVE MODEL FITTING;  SCALAR CASE')                  
!xx    4 FORMAT(  '   < AUTOREGRESSIVE MODEL >',/,1H ,10X,'Z(I) = A(1)*Z(I-
!xx     11) + A(2)*Z(I-2) +  ...  + A(M)*Z(I-M) + E(I)',/,'   WHERE',/,11X,
!xx     2'M:     ORDER OF THE MODEL',/,11X,'E(I):  GAUSSIAN WHITE NOISE WIT
!xx     3H MEAN 0  AND  VARIANCE SD(M).' )                                 
end
