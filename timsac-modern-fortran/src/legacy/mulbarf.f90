! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine mulbarf( zs,n,id,c,lag,zmean,zvari,sd,aic,dic,imin,&
&aicm,sdmin,bw1,bw2,a,b,g,h,e,aicb )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM  MULBAR                                                   
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
!     TIMSAC 78.2.2.                                                    
!     ___          _                  __                                
!     MULTIVARIATE BAYESIAN METHOD OF AR MODEL FITTING                  
!                                                                       
!     THIS PROGRAM DETERMINES MULTI-VARIATE AUTOREGRESSIVE MODELS BY A  
!     BAYESIAN PROCEDURE.  THE BASIC LEAST SQUARES ESTIMATES OF THE PARA
!     ARE OBTAINED BY THE HOUSEHOLDER TRANSFORMATION.                   
!                                                                       
!     THE STATISTIC AIC IS DEFINED BY                                   
!                                                                       
!            AIC  =  N * LOG( DET(SD) ) + 2 * (NUMBER OF PARAMETERS)    
!                                                                       
!       WHERE                                                           
!           N:    NUMBER OF DATA,                                       
!           SD:   ESTIMATE OF INNOVATION VARIANCE MATRIX                
!           DET:  DETERMINANT,                                          
!           K:    NUMBER OF FREE PARAMETERS.                            
!                                                                       
!     BAYESIAN WEIGHT OF THE M-TH ORDER MODEL IS DEFINED BY             
!         W(M)  = CONST * C(M) / (M+1)                                  
!     WHERE                                                             
!         CONST = NORMALIZING CONSTANT                                  
!         C(M)  = EXP( -0.5*AIC(M) ).                                   
!     THE BAYESIAN ESTIMATES OF PARTIAL AUTOREGRESSION COEFFICIENT MATRI
!     OF FORWARD AND BACKWARD MODELS ARE OBTAINED BY (M=1,...,LAG)      
!         G(M)  = G(M)*D(M)                                             
!         H(M)  = H(M)*D(M),                                            
!     WHERE THE ORIGINAL G(M) AND H(M) ARE THE (CONDITIONAL) MAXIMUM    
!     LIKELIHOOD ESTIMATES OF THE HIGHEST ORDER COEFFICIENT MATRICES OF 
!     FORWARD AND BACKWARD AR MODELS OF ORDER M AND D(M) IS DEFINED BY  
!         D(M)  = W(M) + ... + W(LAG).                                  
!                                                                       
!     THE EQUIVALENT NUMBER OF PARAMETERS FOR THE BAYESIAN MODEL IS     
!     DEFINED BY                                                        
!         EK = (D(1)**2 + ... + D(LAG)**2)*ID + ID*(ID+1)/2             
!     WHERE ID DENOTES DIMENSION OF THE PROCESS.                        
!                                                                       
!                                                                       
!       --------------------------------------------------------------- 
!       REFERENCES:                                                     
!          H.AKAIKE(1978), "A BAYESIAN EXTENSION OF THE MINIMUM AIC     
!          PROCEDURE OF AUTOREGRESSIVE MODEL FITTING.",  RESEARCH MEMO. 
!          NO. 126, THE INSTITUTE OF STATISTICAL MATHEMATICS; TOKYO.    
!                                                                       
!          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
!          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
!          30,B,351-363.                                                
!       --------------------------------------------------------------- 
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
!             MRDATA                                                    
!             MREDCT                                                    
!             MARFIT                                                    
!       ----------------------------------------------------------------
!       INPUTS REQUIRED:                                                
!           MT:    INPUT DEVICE FOR ORIGINAL DATA (MT=5; CARD READER)   
!           LAG:   UPPER LIMIT OF AR-ORDER,  MUST BE LESS THAN 31       
!                                                                       
!-----  THE FOLLOWING INPUTS ARE REQUIRED AT SUBROUTINE MRDATA  -----   
!                                                                       
!           TITLE: SPECIFICATION OF DATA                                
!           N:     DATA LENGTH,  MUST BE LESS THAN OR EQUAL TO 1000     
!           ID:    DIMENSION OF VECTOR,  MUST BE LESS THAN 11           
!                       < ID * LAG  MUST BE LESS THAN 101 >             
!           IFM:   CONTROL FOR INPUT                                    
!           FORM:  INPUT DATA FORMAT SPECIFICATION STATEMENT            
!                  -- FOR EXAMPLE --     (8F10.5)                       
!           C(I):  CALIBRATION OF CHANNEL I (I=1,ID)                    
!           Z:     ORIGINAL DATA; Z(K,I) (K=1,N) REPRESENTS THE I-TH CHA
!                  RECORD                                               
!
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4  Z                                                       
!c      DIMENSION  Z(1500,10)                                             
!c      DIMENSION  X(200,100) , D(200)                                    
!c      DIMENSION  A(10,10,30) , B(10,10,30) , G(10,10,30) , H(10,10,30)  
!c      DIMENSION  E(10,10)                                               
!xx      DIMENSION  Z(N,ID), ZS(N,ID), C(ID)
!xx      DIMENSION  ZMEAN(ID), ZVARI(ID)
!xx      DIMENSION  X((LAG+1)*ID*2,(LAG+1)*ID)
!xx      DIMENSION  A(ID,ID,LAG), B(ID,ID,LAG), G(ID,ID,LAG), H(ID,ID,LAG)
!xx      DIMENSION  E(ID,ID)
!
!xx      DIMENSION  SD(LAG+1), AIC(LAG+1), DIC(LAG+1)
!xx      DIMENSION  BW1(LAG+1), BW2(LAG)
integer n, id, lag, imin
real(dp) zs(n,id), c(id), zmean(id), zvari(id), sd(lag+1),&
&aic(lag+1), dic(lag+1), aicm, sdmin, bw1(lag+1),&
&bw2(lag), a(id,id,lag), b(id,id,lag),&
&g(id,id,lag), h(id,id,lag), e(id,id), aicb
! local
integer ipr, ksw, mj, mj1, mj2, n0, nmk
real(dp) z(n,id), x((lag+1)*id*2,(lag+1)*id), ek
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
!c      MJ = 1500                                                         
!c      MJ1 = 200                                                         
!c      MJ2 = 10                                                          
mj = n
mj1 = (lag+1)*id*2
mj2 = id
ipr = 1
ipr = 3
ipr = 2
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      READ( 5,1 )     LAG                                               
!c      WRITE( 6,3 )                                                      
!c      WRITE( 6,4 )                                                      
!c      WRITE( 6,5 )     LAG , MT                                         
!                                                                       
!     --  ORIGINAL DATA LOADING AND MEANS DELETION  --                  
!                                                                       
!c      CALL  MRDATA( MT,MJ,Z,N,ID )                                      
call mrdata( zs,z,n,id,c,zmean,zvari )
!c      CLOSE( MT )
n0 = 0
nmk = n - lag
ksw = 0
!                                                                       
!     --  HOUSEHOLDER REDUCTION  --                                     
!                                                                       
!c      CALL  MREDCT( Z,D,NMK,N0,LAG,ID,MJ,MJ1,KSW,X )                    
x(1:mj1,1:(lag+1)*id) = 0.0d0
call  mredct( z,nmk,n0,lag,id,mj,mj1,ksw,x )
!                                                                       
!     --  AR-MODEL FITTING (BAYESIAN PROCEDURE)  --                     
!                                                                       
!c      CALL  MBYSAR( X,D,NMK,LAG,ID,KSW,IPR,MJ1,MJ2,A,B,G,H,E,AIC,EK )   
!xx      CALL  MBYSAR( X,NMK,LAG,ID,KSW,IPR,MJ1,MJ2,SD,AIC,DIC,
call  mbysar( x,nmk,lag,id,ksw,mj1,mj2,sd,aic,dic,&
&aicm,sdmin,imin,bw1,bw2,a,b,g,h,e,aicb,ek )
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
return
!xx    1 FORMAT( 16I5 )                                                    
!xx    3 FORMAT( 1H ,'PROGRAM TIMSAC 78.2.2',/,'   EXPONENTIALLY WEIGHTED B
!xx     1AYESIAN AUTOREGRESSIVE MODEL FITTING;  MULTI-VARIATE CASE' )      
!xx    4 FORMAT( 1H ,'  < AUTOREGRESSIVE MODEL >',/,1H ,10X,'Z(I) = A(1)*Z(
!xx     1I-1) + A(2)*Z(I-2) + ... + A(II)*Z(I-II) + ... + A(M)*Z(I-M) + W(I
!xx     2)',/,'   WHERE',/,11X,'M:     ORDER OF THE MODEL',/,11X,'W(I):  ID
!xx     3-DIMENSIONAL GAUSSIAN WHITE NOISE WITH MEAN 0 AND VARIANCE',
!xx     4' MATRIX E(M).' )                                                          
!xx    5 FORMAT( 1H ,I4,'-TH ORDER BAYESIAN MODEL IS FITTED',/,1H ,2X,'ORIG
!xx     1INAL DATA INPUT DEVICE  MT =',I4 )                                
end
