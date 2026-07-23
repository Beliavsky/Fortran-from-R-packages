! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine  mulmarf( zs,n,id,c,lag,zmean,zvari,sd1,aic1,dic1,im,&
&aicm,sdm,npr,jndf,af,ex,aic,ei,bi,e,b,lmax,aics )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM  MULMAR                                                   
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
!     TIMSAC 78.2.1.                                                    
!     ___                  _                     __                     
!     MULTIVARIATE CASE OF MINIMUM AIC METHOD OF AR MODEL FITTING.      
!                                                                       
!     THIS PROGRAM FITS A MULTI-VARIATE AUTOREGRESSIVE MODEL BY THE MINI
!     AIC PROCEDURE.  ONLY THE POSSIBILITIES OF ZERO COEFFICIENTS AT THE
!     BEGINNING AND END OF THE MODEL ARE CONSIDERED. THE LEAST SQUARES E
!     OF THE PARAMETERS ARE OBTAINED BY THE HOUSEHOLDER TRANSFORMATION. 
!     AIC IS DEFINED BY                                                 
!                                                                       
!            AIC  =  N * LOG( DET(SD) ) + 2 * (NUMBER OF PARAMETERS)    
!                                                                       
!       WHERE                                                           
!           N:    NUMBER OF DATA,                                       
!           SD:   ESTIMATE OF INNOVATION VARIANCE MATRIX                
!           DET:  DETERMINANT,                                          
!           K:    NUMBER OF FREE PARAMETERS.                            
!                                                                       
!                                                                       
!       --------------------------------------------------------------- 
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM.  
!           MRDATA                                                      
!           MREDCT                                                      
!           MARFIT                                                      
!       --------------------------------------------------------------- 
!       REFERENCE:                                                      
!          G.KITAGAWA AND H.AKAIKE(1978), "A PROCEDURE FOR THE MODELING 
!          OF NON-STATIONARY TIME SERIES.",  ANN. INST. STATIST. MATH., 
!          30,B,351-363.                                                
!       --------------------------------------------------------------- 
!       INPUTS REQUIRED:                                                
!           MT:    INPUT DEVICE FOR ORIGINAL DATA (MT=5; CARD READER)   
!           LAG:   UPPER LIMIT OF AR-ORDER,  MUST BE LESS THAN 31       
!                                                                       
!          .....  FOLLOWING INPUTS ARE REQUIRED AT SUBROUTINE MREDCT  ..
!           TITLE: SPECIFICATION OF DATA                                
!           N:     DATA LENGTH,  MUST BE LESS THAN OR EQUAL TO 1000     
!           ID:    DIMENSION OF DATA,  MUST BE LESS THAN 11             
!                       < ID*(M+1) MUST BE LESS THAN 101 >              
!           IFM:   CONTROL FOR INPUT                                    
!           FORM:  INPUT DATA FORMAT SPECIFICATION STATEMENT            
!                  -- FOR EXAMPLE --     (8F10.5)                       
!          C(I):  CALIBRATION OF CHANNEL I (I=1,ID)                     
!           Z:    ORIGINAL DATA; Z(K,I) (K=1,N) REPRESENTS THE I-TH CHAN
!                 RECORD                                                
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4  Z                                                       
!c      DIMENSION  Z(1500,5)                                              
!c      DIMENSION  X(200,100) , D(200)                                    
!c      DIMENSION  Y(100,100) , B(10,10,30) , E(10,10)                    
!c      DIMENSION  C(10) , EX(10)                                         
!xx      DIMENSION  Z(N,ID), ZS(N,ID), C(ID)
!xx      DIMENSION  ZMEAN(ID), ZVARI(ID)
!xx      DIMENSION  X((LAG+1)*ID*2,(LAG+1)*ID)
!xx      DIMENSION  B(ID,ID,LAG) , E(ID,ID), BI(ID,ID,LAG) , EI(ID,ID)
!xx      DIMENSION  EX(ID), CV(ID)
!xx      DIMENSION  SD1(LAG+1,ID), AIC1(LAG+1,ID), DIC1(LAG+1,ID)
!xx      DIMENSION  AICM(ID), SDM(ID), IM(ID)
!xx      DIMENSION  JNDF((LAG+1)*ID,ID), AF((LAG+1)*ID,ID)
!xx      DIMENSION  NPR(ID), AIC(ID)
integer n, id, lag, im(id), npr(id), jndf((lag+1)*id,id), lmax
real(dp) zs(n,id), c(id), zmean(id), zvari(id),&
&sd1(lag+1,id), aic1(lag+1,id), dic1(lag+1,id),&
&aicm(id), sdm(id), af((lag+1)*id,id), ex(id),&
&aic(id), ei(id,id), bi(id,id,lag), e(id,id),&
&b(id,id,lag), aics
! local
integer ipr, ksw, mj, mj1, mj2, mj3, mj4, n0, nmk
real(dp) z(n,id), x((lag+1)*id*2,(lag+1)*id), cv(id)
!x      INTEGER*1  TMP(1)
!x      CHARACTER  CNAME*80
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
!           MJ:    ABSOLUTE DIMENSION FOR SUBROUTINE CALL               
!           MJ1:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL, SHOULD BE LAR
!                  THAN ID*(M+1)                                        
!           MJ2:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL               
!           MJ3:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL               
!           MJ4:   ABSOLUTE DIMENSION FOR SUBROUTINE CALL               
!           IPR:   PRINT OUT CONTROL                                    
!                                                                       
!c      MJ = 1500                                                         
!c      MJ1 = 200                                                         
!c      MJ2 = 10                                                          
!c      MJ3 = 30                                                          
!c      MJ4 = 100                                                         
mj = n
mj1 = (lag+1)*id*2
mj2 = id
mj3 = lag
mj4 = (lag+1)*id
ksw = 0
ipr = 3
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      READ( 5,1 )     LAG                                               
!c      WRITE( 6,3 )                                                      
!c      WRITE( 6,4 )     LAG , MT                                         
!                                                                       
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
!xcx            WRITE(*,*) ' ***  mulmar temp FILE OPEN ERROR :',CNAME,IVAR
!x            IER=IVAR
!x            IFG=0
!x         END IF
!x      END IF

!
!     --  ORIGINAL DATA LOADING AND MEANS DELETION  --                  
!                                                                       
!c      CALL  MRDATA( MT,MJ,Z,N,ID )                                      
call mrdata( zs,z,n,id,c,zmean,zvari )
!c      CLOSE( MT )
n0 = 0
nmk = n - lag
!                                                                       
!     --  HOUSEHOLDER REDUCTION  --                                     
!                                                                       
!c      CALL  MREDCT( Z,D,NMK,N0,LAG,ID,MJ,MJ1,KSW,X )
x(1:mj1,1:mj4) = 0.0d0
call  mredct( z,nmk,n0,lag,id,mj,mj1,ksw,x )
!                                                                       
!     --  AR-MODEL FITTING (MAICE PROCEDURE)  --                        
!                                                                       
!c      CALL  MARFIT( X,Y,D,NMK,ID,LAG,KSW,MJ1,MJ2,MJ3,MJ4,0,IPR,B,E,EX,C,
!c     *              LMAX,AIC )                                          
call marfit( x,nmk,id,lag,ksw,mj1,mj2,mj3,mj4,0,ipr,aic1,sd1,dic1,&
!x     * AICM,SDM,IM,BI,EI,B,E,EX,CV,LMAX,AICS,JNDF,AF,NPR,AIC,IFG,LU )
&aicm,sdm,im,bi,ei,b,e,ex,cv,lmax,aics,jndf,af,npr,aic )
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
!c      CLOSE( LU )
!x      IF( IFG.NE.0 ) CLOSE( LU )
return
!xx    1 FORMAT( 16I5 )                                                    
!xx    3 FORMAT( ' PROGRAM TIMSAC 78.2.1',/'   MULTI-VARIATE AUTOREGRESSIVE
!xx     * MODEL FITTING  ;  LEAST SQUARES METHOD BY HOUSEHOLDER TRANSFORMAT
!xx     2ION',/,'   < AUTOREGRESSIVE MODEL >',/,1H ,10X,'Z(I) = A(1)*Z(I-1)
!xx     3 + A(2)*Z(I-2) + ... + A(II)*Z(I-II) + ... + A(M)*Z(I-M) + W(I)',/
!xx     4,'   WHERE',/,11X,'M:     ORDER OF THE MODEL',/,11X,'W(I):  ID-DIM
!xx     5ENSIONAL GAUSSIAN WHITE NOISE WITH MEAN 0 AND VARIANCE MATRIX E(M)
!xx     6.' )                                                              
!xx    4 FORMAT( 1H ,'FITTING UP TO THE ORDER  K =',I3,'   IS TRIED',/,' OR
!xx     1IGINAL DATA INPUT DEVICE   MT =',I3 )                             
end
