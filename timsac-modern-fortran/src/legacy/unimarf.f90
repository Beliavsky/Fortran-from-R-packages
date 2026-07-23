! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine  unimarf( zs,n,lag,zmean,sum,sd,aic,dic,m,aicm,sdm,a )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM  UNIMAR                                                   
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
!     TIMSAC 78.1.1.                                                    
!     ___                _                     __                       
!     UNIVARIATE CASE OF MINIMUM AIC METHOD OF AR MODEL FITTING.        
!                                                                       
!       THIS IS THE BASIC PROGRAM FOR THE FITTING OF AUTOREGRESSIVE MODE
!       OF SUCCESSIVELY HIGHER ORDERS BY THE METHOD OF LEAST SQUARES    
!       REALIZED THROUGH HOUSEHOLDER TRANSFORMATION.  THE OUTPUTS ARE   
!       THE ESTIMATES OF THE COEFFICIENTS, THE INNOVATION VARIANCES AND 
!       CORRESPONDING AIC STATISTICS.  AIC IS DEFINED BY                
!                                                                       
!               AIC  =  N * LOG( SD )  +  2 * ( NUMBER OF PARAMETERS )  
!       WHERE                                                           
!             N:    DATA LENGTH,                                        
!             SD:   ESTIMATE OF THE INNOVATION VARIANCE.                
!                                                                       
!       --------------------------------------------------------------- 
!       REFERENCE:                                                      
!          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
!          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
!          30,B,351-363.                                                
!       --------------------------------------------------------------- 
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
!             REDATA                                                    
!             REDUCT                                                    
!             ARMFIT                                                    
!       --------------------------------------------------------------- 
!       INPUTS REQUIRED:                                                
!             MT:    INPUT DEVICE SPECIFICATION (MT=5 : CARD READER)    
!             LAG:   UPPER LIMIT OF AR-ORDER, MUST BE LESS THAN 101     
!                                                                       
!       --  THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDATA  -- 
!                                                                       
!             TITLE: TITLE OF DATA                                      
!             N:     DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000   
!             DFORM: INPUT DATA FORMAT SPECIFICATION STATEMENT          
!                    -- FOR EXAMPLE --     (8F10.5)                     
!             (Z(I),I=1,N):  ORIGINAL DATA                              
!       --------------------------------------------------------------- 
!
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4   Z(10000) , TITLE(20)                                   
!c      REAL * 4   TITLE(20)
!c      DIMENSION  Z(10000)
!c      DIMENSION  X(200,101) , D(200) , A(100)                           
!xx      DIMENSION  ZS(N), Z(N)
!xx      DIMENSION  X(N+1,LAG+1), A(LAG)
!xx      DIMENSION  SD(LAG+1), AIC(LAG+1), DIC(LAG+1)
integer n, lag, m
real(dp) zs(n), zmean, sum, sd(lag+1), aic(lag+1),&
&dic(lag+1), aicm, sdm, a(lag)
! local
integer isw, k, mj1, nmk
real(dp) z(n), x(n+1,lag+1)
!x      INTEGER*1  TMP(1)
!x      CHARACTER  CNAME*80
!                                                                       
!        EXTERNAL SUBROUTINE DECLARATION:                               
!                                                                       
external  setx1
!
!c      CHARACTER(100) IFLNAM,OFLNAM
!c      CALL FLNAM2( IFLNAM,OFLNAM,NFL )
!c      IF ( NFL.EQ.0 ) GO TO 999
!c      IF ( NFL.EQ.2 ) THEN
!c         OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR )
!c      ELSE
!c         CALL SETWND
!c      END IF
!x      IER=0
!x      LU=3
!x      DO 100 I = 1,80
!x         CNAME(I:I) = ' '
!x  100 CONTINUE
!x      I = 1
!x      IFG = 1
!x      DO WHILE( (IFG.EQ.1) .AND. (I.LE.80) )
!x	   IF ( TMP(I).NE.ICHAR(' ') ) THEN
!x            CNAME(I:I) = CHAR(TMP(I))
!x            I = I+1
!x         ELSE
!x            IFG = 0
!x         END IF
!x      END DO
!x      IF ( I.GT.1 ) THEN
!x         IFG = 1
!x         OPEN (LU,FILE=CNAME,IOSTAT=IVAR)
!x         IF (IVAR .NE. 0) THEN
!xcx            WRITE(*,*) ' ***  unimar temp FILE OPEN ERROR :',CNAME,IVAR
!x            IER=IVAR
!x            IFG=0
!x         END IF
!x      END IF
!                                                                       
!        PARAMETERS:                                                    
!             MJ1:  ABSOLUTE DIMENSION FOR SUBROUTINE CALL              
!             ISW:  =1   PARAMETERS OF MAICE MODEL ONLY ARE REQUESTED   
!                   =2   PARAMETERS OF ALL MODELS ARE REQUESTED         
!                                                                       
!c      MJ1 = 200                                                         
mj1 = n+1
isw = 2
!                                                                       
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      READ( 5,1 )     LAG                                               
!c      WRITE( 6,3 )                                                      
!c      WRITE( 6,2 )     LAG                                              
!c      WRITE( 6,4 )     MT                                               
!                                                                       
!          +-----------------------------------------+                +-
!          ! ORIGINAL DATA LOADING AND MEAN DELETION !                ! 
!          +-----------------------------------------+                +-
!                                                                       
!c      CALL  REDATA( Z,N,MT,TITLE )                                      
!c	CLOSE( MT )
call  redata( zs,z,n,zmean,sum )
nmk = n - lag
k = lag
!                                                                       
!          +-----------------------+                                  +-
!          ! HOUSEHOLDER REDUCTION !                                  ! 
!          +-----------------------+                                  +-
!                                                                       
!c      CALL  REDUCT( SETX1,Z,D,NMK,0,K,MJ1,LAG,X )                       
call  reduct( setx1,z,nmk,0,k,mj1,lag,x )
!                                                                       
!          +------------------+                                       +-
!          ! AR MODEL FITTING !                                       ! 
!          +------------------+                                       +-
!                                                                       
!c      CALL  ARMFIT( X,K,LAG,NMK,ISW,TITLE,MJ1,A,SD,M )                  
!x      CALL  ARMFIT( X,K,LAG,NMK,ISW,MJ1,A,M,SD,AIC,DIC,SDM,AICM,
!x     *              IFG,LU )
call  armfit( x,k,lag,nmk,isw,mj1,a,m,sd,aic,dic,sdm,aicm )
!
!c      GO TO 999
!                                                                      +
!c  900 CONTINUE
!c      WRITE(6,600) IVAR,OFLNAM
!c  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
!c      GO TO 999
!                                                                      !
!c  910 CONTINUE
!c      IF ( NFL.EQ.2 ) CLOSE( 6 )
!c#ifdef __linux__
!cC      reopen #6 as stdout
!c      IF ( NFL.EQ.2 ) OPEN(6, FILE='/dev/fd/1')
!c#endif
!cC /* __linux__ */
!c      WRITE(6,610) IVAR,IFLNAM
!c  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
!                                                                      +
!c  999 CONTINUE
!x      IF (IFG.NE.0) CLOSE(LU)                                           
return
!                                                                       
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( 1H ,'FITTING UP TO THE ORDER  K =',I3,'  IS TRIED' )      
!xx    3 FORMAT( ' PROGRAM TIMSAC 78.1.1',/'   AUTOREGRESSIVE MODEL FITTING
!xx     1  (SCALAR CASE)  ;   LEAST SQUARES METHOD BY HOUSEHOLDER TRANSFORM
!xx     2ATION',/,'   < AUTOREGRESSIVE MODEL >',/,1H ,10X,'Z(I) = A(1)*Z(I-
!xx     31) + A(2)*Z(I-2) +  ...  + A(M)*Z(I-M) + E(I)',/,'   WHERE',/,11X,
!xx     4'M:     ORDER OF THE MODEL',/,11X,'E(I):  GAUSSIAN WHITE NOISE WIT
!xx     5H MEAN 0  AND  VARIANCE SD(M).' )                                 
!xx    4 FORMAT( 1H ,'ORIGINAL DATA INPUT DEVICE  MT =',I3 )               
!                                                                       
end
