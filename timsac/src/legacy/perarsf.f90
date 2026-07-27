! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine perarsf( zs,n,ip,lag,ksw,zmean,sum,npr,jndf,af,aicf,&
&b,e,c,ex,lmax)
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM  PERARS                                                   
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
!     TIMSAC 78.2.3.                                                    
!     ___      _   _                _                                   
!     PERIODIC AUTOREGRESSION FOR A SCALAR TIME SERIES                  
!                                                                       
!       THIS IS THE PROGRAM FOR THE FITTING OF PERIODIC AUTOREGRESSIVE M
!       BY THE METHOD OF LEAST SQUARES REALIZED THROUGH HOUSEHOLDER     
!       TRANSFORMATION.  THE OUTPUTS ARE THE ESTIMATES OF THE REGRESSION
!       COEFFICIENTS AND INNOVATION VARIANCE OF THE PERIODIC AR-MODEL FO
!       EACH INSTANT.                                                   
!                                                                       
!       PERIOD I           1              2                             
!                    !-----------!  !-----------!                  !----
!            Z(II)   +--+--+--+--+--+--+--+--+--+   . . . . . .    +--+-
!            Y(I,J)  +--+--+--+--+  +--+--+--+--+                  +--+-
!       INSTANT J    1  2       IP  1  2       IP                  1  2 
!                                                                       
!       WHERE                                                           
!          IP:     NUMBER OF INSTANTS IN ONE PERIOD                     
!          ND:     NUMBER OF PERIODS                                    
!          Y(I,J) = Z(IP*(I-1)+J)                                       
!                                                                       
!       THE STATISTIC AIC IS DEFINED BY                                 
!                                                                       
!               AIC  =  N * LOG( SD )  +  2 * ( NUMBER OF PARAMETERS )  
!       WHERE                                                           
!             N:    DATA LENGTH,                                        
!             SD:   ESTIMATE OF THE INNOVATION VARIANCE.                
!                                                                       
!                                                                       
!       --------------------------------------------------------------- 
!       REFERRENCES:                                                    
!          R.H.JONES AND W.M.BRELSFORD(1967), "TIME SERIES WITH PERIODIC
!          STRUCTURE.",  BIOMETRIKA,54,403-408.                         
!                                                                       
!          M.PAGANO(1978), "ON PERIODIC AND MULTIPLE AUTOREGRESSIONS."  
!          ANN. STATIST., 6, 1310-1317.                                 
!       --------------------------------------------------------------- 
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM.  
!            REDATA                                                     
!            PERREG                                                     
!            MREDCT                                                     
!            MARFIT                                                     
!            PRINT4                                                     
!       --------------------------------------------------------------- 
!       INPUTS REQUIRED:                                                
!             MT:    INPUT DEVICE SPECIFICATION (MT=5 : CARD READER)    
!             IP:    NUMBER OF OBSERVATIONS WITHIN A PERIOD             
!             LAG:   MAXIMUM LAG OF PERIODS                             
!                                                                       
!       --  THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDATA  -- 
!                                                                       
!             TITLE: TITLE OF DATA                                      
!             N:     DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000   
!             DFORM: INPUT DATA FORMAT SPECIFICATION STATEMENT          
!                    -- EXAMPLE --    (8F10.5)                          
!             (Z(I),I=1,N):  ORIGINAL DATA                              
!       --------------------------------------------------------------- 
!
!xx      IMPLICIT  REAL * 8 ( A-H,O-Z )                                    
!C      REAL  * 4   Z , Y                                                 
!c      REAL * 4   TITLE(20)                                              
!c      DIMENSION  Z(5000) , Y(200,24)                                    
!c      DIMENSION  X(200,150) , D(300)
!c      DIMENSION  U(150,150) , B(24,24,5) , E(24,24)                     
!c      DIMENSION  C(24)                                                  
!c      DIMENSION  EX(24)                                                 
!xx      DIMENSION  ZS(N), Z(N), Y(N/IP,IP)                                    
!xx      DIMENSION  X(((LAG+1)*IP+KSW)*2,(LAG+1)*IP+KSW)
!xx      DIMENSION  B(IP,IP,LAG) , E(IP,IP), BI(IP,IP,LAG) , EI(IP,IP)                     
!xx      DIMENSION  C(IP)
!xx      DIMENSION  EX(IP)
!xx      DIMENSION  AIC(LAG+1,IP), SD(LAG+1,IP), DIC(LAG+1,IP)
!xx      DIMENSION  AICM(IP), SDM(IP), IM(IP)
!xx      DIMENSION  JNDF((LAG+1)*IP+KSW,IP), AF((LAG+1)*IP+KSW,IP)
!xx      DIMENSION  NPR(IP), AICF(IP)
!
integer n, ip, lag, ksw, npr(ip), jndf((lag+1)*ip+ksw,ip), lmax
real(dp) zs(n), zmean, sum, af((lag+1)*ip+ksw,ip),&
&aicf(ip), b(ip,ip,lag), e(ip,ip), c(ip), ex(ip)
! local
integer id, ipr, isw, mj, mj1, mj2, mj3, mj4, n0, nd, nmk, im(ip)
real(dp) x(((lag+1)*ip+ksw)*2,(lag+1)*ip+ksw), y(n/ip,ip),&
&z(n), bi(ip,ip,lag), ei(ip,ip), aic(lag+1,ip),&
&sd(lag+1,ip), dic(lag+1,ip), aicm(ip), sdm(ip),&
&aics
!
!c      CHARACTER(100) IFLNAM,OFLNAM
!c      CALL FLNAM2( IFLNAM,OFLNAM,NFL )
!c      IF ( NFL.EQ.0 ) GO TO 999
!c      IF ( NFL.EQ.2 ) THEN
!c         OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR )
!c      ELSE
!c         CALL SETWND
!c      END IF
!                                                                       
!
!          PARAMETERS:                                                  
!                                                                       
!c      MJ = 200                                                          
!c      MJ1 = 200                                                         
!c      MJ2 = 24                                                          
!c      MJ3 = 5                                                           
!c      MJ4 = 150                                                         
mj = n/ip
mj2 = ip
mj3 = lag
mj4 = (lag+1)*ip+ksw
mj1 = mj4*2
isw = 1
ipr = 2
!                                                                       
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      READ( 5,1 )     IP , LAG , KSW                                    
!c      WRITE( 6,2 )                                                      
!c      WRITE( 6,3 )                                                      
!c      WRITE( 6,4 )   IP , LAG , KSW , MT                                
!c      WRITE( 6,5 )                                                      
!                                                                       
!          ORIGINAL DATA LOADING                                        
!                                                                       
!c      CALL  REDATA( Z,N,MT,TITLE )                                      
call  redata( zs,z,n,zmean,sum )
!c      CLOSE( MT )
!                                                                       
!          DATA MATRIX SET UP                                           
!                                                                       
call  perreg( z,n,ip,mj,y,nd )
nmk = nd - lag
n0 = 0
id = ip
!                                                                       
!          REDUCTION TO AN UPPER TRIANGULAR FORM                        
!                                                                       
!c      CALL  MREDCT( Y,D,NMK,N0,LAG,ID,MJ,MJ1,KSW,X )
x(1:mj1,1:mj4) = 0.0d0
call  mredct( y,nmk,n0,lag,id,mj,mj1,ksw,x )
!                                                                       
!          INSTANTANEOUS RESPONSE MODEL FITTING                         
!                                                                       
!c      CALL  MARFIT( X,U,D,NMK,IP,LAG,KSW,MJ1,MJ2,MJ3,MJ4,ISW,IPR,B,E,EX,
!c     *              C,LMAX,AIC )                                        
!x      IFG = 0
call marfit( x,nmk,ip,lag,ksw,mj1,mj2,mj3,mj4,isw,ipr,aic,sd,&
!x     *DIC,AICM,SDM,IM,BI,EI,B,E,EX,C,LMAX,AICS,JNDF,AF,NPR,AICF,IFG,LU )
&dic,aicm,sdm,im,bi,ei,b,e,ex,c,lmax,aics,jndf,af,npr,aicf )
!   
!                                                                       
!          REGRESSION MODEL PRINT OUT                                   
!                                                                       
!c      CALL  PRINT4( B,E,C,EX,ID,LMAX,MJ2 )                              
!c      GO TO 999
!                                                                       
!c  900 CONTINUE
!c      WRITE(6,600) IVAR,OFLNAM
!c  600 FORMAT(/,' !!! Output_Data_File OPEN ERROR ',I8,/,5X,100A)
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
!c  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,/,5X,100A)
!
!c  999 CONTINUE
return
!                                                                       
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( 1H ,'PROGRAM  TIMSAC 78.2.3.',/,1H ,'  PERIODIC AUTOREGRES
!xx     1SIVE MODELS FITTING BY THE METHOD OF LEAST SQUARES (SCALAR CASE)')
!xx    3 FORMAT( 1H ,'  <PERIODIC AUTOREGRESSIVE MODEL (J=1,...,IP) >',/,  
!xx     1 1H ,8X,'Y(I,J) = C(J) + A(1,J,0)*Y(I,1) + ... + A(J-1,J,0)*',
!xx     2'Y(I,J-1) + A(1,J,1)*Y(I-1,1) + ... + A(IP,J,1)*Y(I-1,IP)',
!xx     3' + ... + E(I,J)',/,   1H ,'  WHERE',/,
!xx     4  1H ,7X,'IP:       NUMBER OF INSTANTS IN ONE PERIOD',/,          
!xx     5  1H ,7X,'E(I,J):   GAUSSIAN WHITE NOISE' )                       
!xx    4 FORMAT( 1H ,2X,'IP =',I3,5X,'LAG =',I3,5X,'KSW =',I3,/,3X,'ORIGINA
!xx     *L DATA INPUT DEVICE   MT =',I3 )                                  
!xx    5 FORMAT( 1H ,'  *****  WHEN KSW IS SET TO 0, THE CONSTANT TERM C(J)
!xx     * IS EXCLUDED.  *****' )                                           
end
subroutine  perreg( z,n,ip,mj,y,nd )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PREPARES DATA MATRIX Y FOR THE FITTING OF A PERIOD
!     AUTOREGRESSIVE MODEL FROM THE DATA VECTOR Z OF CONSECUTIVE OBSERVA
!     EACH COLUMN OF Y IS COMPOSED OF THE OBSERVATIONS AT THE SAME INSTA
!     WITHIN A PERIOD.                                                  
!                                                                       
!       INPUTS:                                                         
!          Z:     ORIGINAL DATA VECTOR                                  
!          N:     LENGTH OF ORIGINAL DATA                               
!          IP:    SPAN OF ONE PERIOD                                    
!          MJ:    ABSOLUTE DIMENSION OF Y                               
!                                                                       
!       OUTPUTS:                                                        
!          Y:     REARRANGED DATA MATRIX                                
!          ND:    NUMBER OF ROWS OF Y                                   
!                                                                        
!C      DIMENSION  Z(1) , Y(MJ,1)                                         
!x      REAL * 8  Z(1) , Y(MJ,1)
!xx      REAL * 8  Z(N) , Y(MJ,IP)
integer n, ip, mj, nd
real(dp) z(n), y(mj,ip)
! local
integer i, ii, j
nd = n / ip
!xx      DO 10  I=1,ND
do 20  i=1,nd
do 10  j=1,ip
ii = (i-1)*ip + j
!xx   10 Y(I,J) = Z(II)
y(i,j) = z(ii)
10 continue
20 continue
return
!                                                                       
end
