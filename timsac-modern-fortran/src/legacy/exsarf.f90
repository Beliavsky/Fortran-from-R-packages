! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine  exsarf( z1,n,lag,zmean,sum,sd,aic,dic,m1,amin,sdm1,a1,&
!x     *                    SDM2,A2,TMP,IER,JER )
&sdm2,a2,jer )
  use timsac_kinds, only: dp
  implicit none
!----------------   for  isw = 2   --------------------
!xc     *                    SDM22,A22,TMP )
!     *                    SDM22,A22,JER )
!
!c      PROGRAM  EXSAR                                                    
!.......................................................................
!.....PLANNED BY H.AKAIKE...............................................
!.....DESIGNED BY H.AKAIKE AND G.KITAGAWA...............................
!.....PROGRAMMED BY G.KITAGAWA AND F.TADA...............................
!.....ADDRESS: THE INSTITUTE OF STATISTICAL MATHEMATICS, 4-6-7 MINAMI-AZ
!..............MINATO-KU, TOKYO 106, JAPAN..............................
!.....DATE OF THE LATEST REVISION:  JUN. 29, 1979.......................
!.......................................................................
!.....THIS PROGRAM WAS ORIGINALLY PUBLISHED IN "TIMSAC-78", BY H.AKAIKE,
!.....G.KITAGAWA, E.ARAHATA AND F.TADA, COMPUTER SCIENCE MONOGRAPHS, NO.
!.....THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO, 1979.............
!.......................................................................
!     TIMSAC 78.5.1                                                     
!     __                                 _      __                      
!     EXACT MAXIMUM LIKELIHOOD METHOD OF SCALAR AR-MODEL FITTING        
!                                                                       
!     THIS PROGRAM PRODUCES EXACT MAXIMUM LIKELIHOOD ESTIMATES OF THE   
!     PARAMETERS OF A SCALAR AR-MODEL.                                  
!                                                                       
!     THE AR-MODEL IS GIVEN BY                                          
!                                                                       
!               Z(I) = A(1)*Z(I-1) + ... + A(K)*Z(I-K) + E(I)           
!                                                                       
!     WHERE E(I) IS A ZERO MEAN WHITE NOISE.                            
!                                                                       
!     --------------------------------------------------------------    
!     THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:    
!             REDATA                                                    
!             REDUCT                                                    
!             ARMFIT                                                    
!             RECOEF                                                    
!             ARMLE                                                     
!             PRINTA                                                    
!     --------------------------------------------------------------    
!     INPUTS REQUIRED:                                                  
!          MT:      INPUT DEVICE SPECIFICATION (MT=5: CARD READER)      
!          LAG:     UPPER LIMIT OF AR-ORDER, MUST BE LESS THAN 51       
!                                                                       
!     --  THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE REDATA  --   
!          TITLE:  TITLE OF DATA                                        
!          N:      DATA LENGTH, MUST BE LESS THAN OR EQUAL TO 10000     
!          DFORM:  INPUT DATA FORMAT SPECIFICATION STATEMENT            
!                  -- EXAMPLE --     (8F10.5)                           
!          (Z(I),I=1,N):  ORIGINAL DATA                                 
!     ---------------------------------------------------------------   
!                                                                       
!xx      IMPLICIT  REAL *8  ( A-H , O-Z )                                  
!C      REAL * 4  Z(10000) , TITLE(20)                                    
!c      REAL * 4   TITLE(20)
!c      DIMENSION  Z(10000)
!c      DIMENSION  X(200,51) , D(200) , A(50)                             
!xx      DIMENSION  Z1(N), Z(N)
!xx      DIMENSION  X(N-LAG,LAG+1), A1(LAG), A2(LAG)
!----------------    isw = 2   -----------------
!xx      DIMENSION  A22(LAG,LAG), SDM22(LAG)
!-----------------------------------------------
integer n, lag, m1, jer
real(dp) z1(n), zmean, sum, sd(lag+1), aic(lag+1),&
&dic(lag+1), amin, sdm1, a1(lag), sdm2, a2(lag)
! local
integer i, ipr, isw, k, m, mj1, nmk
real(dp) z(n), x(n-lag,lag+1)
!----------------    isw = 2   -----------------
real(dp) a22(lag,lag), sdm22(lag)
!-----------------------------------------------

!C      COMMON     / AAA /  N , Z                                         
!c      COMMON     / AAA /  N
!c      COMMON     / BBB /  Z
!c      COMMON     / CCC /  ISW , IPR                                     
!
!xx      DIMENSION  SD(LAG+1), AIC(LAG+1), DIC(LAG+1)
!x      INTEGER*1  TMP(1)
!x      CHARACTER  CNAME*80
!
!c      CHARACTER(100) IFLNAM,OFLNAM
!                                                                       
!       EXTERNAL SUBROUTINE DECLARATION:                                
!                                                                       
external   setx1
!
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
!xcx            WRITE(*,*) ' ***  exsar temp FILE OPEN ERROR :',CNAME,IVAR
!x            IER=IVAR
!x            IFG=0
!x         END IF
!x      END IF
!                                                                       
!          PARAMETERS:                                                  
!             MJ1:  ABSOLUTE DIMENSION FOR SUBROUTINE CALL              
!             ISW:  =1  PARAMETERS OF MAICE MODEL ONLY ARE REQUESTED    
!                   =2  PARAMETERS OF ALL MODELS ARE REQUESTED          
!             IPR:  PRINT OUT CONTROL IN THE NON-LINEAR OPTIMIZATION PRO
!                                                                       
!c      MJ1 = 200                                                         
mj1 = n-lag
isw = 1
!                                                                       
!C      READ( 5,1 )     MT                                                
!c      MT = 5
!c      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD')
!c      READ( 5,1 )     LAG                                               
!c      WRITE( 6,2 )                                                      
!c      WRITE( 6,3 )   LAG , MT                                           
!                                                                       
!          +-----------------------------------------+                  
!          ! ORIGINAL DATA LOADING AND MEAN DELETION !                  
!          +-----------------------------------------+                  
!                                                                       
!c      CALL  REDATA( Z,N,MT,TITLE )                                      
!c      CLOSE( MT )
call  redata( z1,z,n,zmean,sum )
k = lag
nmk = n - k
!                                                                       
!          +-----------------------+                                    
!          ! HOUSEHOLDER REDUCTION !                                    
!          +-----------------------+                                    
!                                                                       
!c      CALL  REDUCT( SETX1,Z,D,NMK,0,K,MJ1,LAG,X )                       
call  reduct( setx1,z,nmk,0,k,mj1,lag,x )
!                                                                       
!          +------------------+                                         
!          ! AR MODEL FITTING !                                         
!          +------------------+                                         
!                                                                       
!c      CALL  ARMFIT( X,K,LAG,NMK,ISW,TITLE,MJ1,A,SD,M )                  
!x      CALL  ARMFIT( X,K,LAG,NMK,ISW,MJ1,A1,M1,SD,AIC,DIC,SDM1,AMIN,
!x     *              IFG,LU )
call  armfit( x,k,lag,nmk,isw,mj1,a1,m1,sd,aic,dic,sdm1,amin )
do 5  i = 1,k
!xx    5 A2(I) = A1(I)
a2(i) = a1(i)
5 continue
!                                                                       
!          +-------------------------------------------+                
!          ! MAXIMIZATION OF EXACT LIKELIHOOD FUNCTION !                
!          +-------------------------------------------+                
!                                                                       
jer = 0
if( isw .eq. 2 )  go to 10
!                                                                       
ipr = 7
!c      CALL  ARMLE( Z,N,M,K,TITLE,A )                                    
!cx      CALL  ARMLE( Z,N,M1,K,A2,SDM2,ISW,IPR,IFG,LU )
!x      CALL  ARMLE( Z,N,M1,K,A2,SDM2,ISW,IPR,IFG,LU,JER )
!xx      CALL  ARMLE( Z,N,M1,K,A2,SDM2,ISW,IPR,JER )
call  armle( z,n,m1,k,a2,sdm2,isw,jer )
!x      IF( JER .NE. 0 ) RETURN
!x      IF (IFG.NE.0) CLOSE(LU)
return
!                                                                       
10 do 20  m=1,k
!c      CALL  RECOEF( X,M,K,MJ1,A )                                       
call  recoef( x,m,k,mj1,a2 )
ipr = 5
!c   20 CALL  ARMLE( Z,N,M,K,TITLE,A )                                    
!cx      CALL  ARMLE( Z,N,M,K,A2,SDM2,ISW,IPR,IFG,LU )
!x      CALL  ARMLE( Z,N,M,K,A2,SDM2,ISW,IPR,IFG,LU,JER )
!xx      CALL  ARMLE( Z,N,M,K,A2,SDM2,ISW,IPR,JER )
call  armle( z,n,m,k,a2,sdm2,isw,jer )
if( jer .ne. 0 ) return
do 15  i = 1,m
!xx   15 A22(I,M) = A2(I)
a22(i,m) = a2(i)
15 continue
sdm22(m) = sdm2
20 continue
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
!cC     reopen #6 as stdout
!c      IF ( NFL.EQ.2 ) OPEN(6, FILE='/dev/fd/1')
!c#endif
!cC /* __linux__ */
!c      WRITE(6,610) IVAR,IFLNAM
!c  610 FORMAT(/,' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
!                                                                       
!xx    1 FORMAT( 16I5 )                                                    
!xx    2 FORMAT( 1H ,'PROGRAM 78.5.1',/,'   EXACT MAXIMUM LIKELIHOOD METHOD
!xx     * OF AUTOREGRESSIVE MODEL FITTING;  SCALAR CASE',/,'   < BASIC AUTO
!xx     *REGRESSIVE MODEL >',/,1H ,10X, 'Z(I) = A(1)*Z(I-1) + ... + A(K)*Z(
!xx     *I-K) + E(I)',/,1H ,'WHERE',/,11X,'K:     ORDER OF THE MODEL',/,11X
!xx     *,'E(I):  ZERO MEAN WHITE NOISE' )                                 
!xx    3 FORMAT( ///1H ,'  FITTING UP TO THE ORDER  K =',I3,'  IS TRIED',/,
!xx     *1H ,'  ORIGINAL DATA INPUT DEVICE  MT =',I3 )                     
!                                                                       
!c  999 CONTINUE
!x      IF (IFG.NE.0) CLOSE(LU)                                                                       
return
end
!c      SUBROUTINE  ARMLE( Z,N,K,L,TITLE,A )                              
!cx      SUBROUTINE  ARMLE( Z,N,K,L,A,SDM,ISW,IPR,IFG,LU )
!x      SUBROUTINE  ARMLE( Z,N,K,L,A,SDM,ISW,IPR,IFG,LU,JER )
!xx      SUBROUTINE  ARMLE( Z,N,K,L,A,SDM,ISW,IPR,JER )
subroutine  armle( z,n,k,l,a,sdm,isw,jer )
  use timsac_kinds, only: dp
  implicit none
!.....DATE OF THE LATEST REVISION:  JUN. 29, 1979.......................
!                                                                       
!     THIS SUBROUTINE PRODUCES EXACT MAXIMUM LIKELIHOOD ESTIMATES OF THE
!     PARAMETERS OF AN AUTOREGRESSIVE MODEL                             
!                                                                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             DAVIDN                                                    
!             PRINTA                                                    
!       ----------------------------------------------------------------
!          INPUTS:                                                      
!             Z:       ORIGINAL DATA                                    
!             N:       DATA LENGTH                                      
!             K:       ORDER OF THE AR-MODEL                            
!             L:       HIGHEST POSSIBLE ORDER OF THE MODEL              
!             TITLE:   TITLE OF DATA                                    
!                                                                       
!          OUTPUTS:                                                     
!             A:       MAXIMUM LIKELIHOOD ESTIMATES OF AR-COEFFICIENTS  
!                                                                       
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4   Z , TITLE                                              
!c      REAL * 4   TITLE
!c      DIMENSION  Z(1) , A(1) , TITLE(1)                                 
!c      DIMENSION  R(31,31) , C(31)                                       
!xx      DIMENSION  Z(N), A(K)
!xx      DIMENSION  R(K+1,K+1), C(L+1)                                       
integer n, k, l, isw, jer
real(dp) z(n), a(k), sdm
! local
integer i, ihes, ii, imj, j, jj, jmi, k1, ki, kk, kmi, l1, n1,&
&nmi, nml, nmlp1
real(dp) r(k+1,k+1), c(l+1), sum, f0, aic, sd
!C      DIMENSION  TTL(4)                                                 
!c      REAL * 4 TTL(8)
!c      COMMON     / DDD /  R , F , AIC , SD                              
!C      DATA         TTL / 8H  MAXIMU,8HM LIKELI,8HHOOD EST,8HIMATE  . /  
!c      DATA         TTL / 4H  MA,4HXIMU,4HM LI,4HKELI,
!c     1					4HHOOD,4H EST,4HIMAT,4HE  . / 
!                                                                       
!       EXTERNAL FUNCTION DECRALATION                                   
!                                                                       
external  funct
external  hesian
!                                                                       
ihes = 1
!                                                                       
!c      WRITE( 6,3 )                                                      
!x      IF( IFG.NE.0 ) WRITE( LU,3 )  K
nmlp1 = n - l + 1
k1 = k + 1
l1 = l + 1
n1 = n + 1
nml = n - l
do  20  ii=1,l1
j = ii - 1
sum = 0.d0
do 10  i=l1,nml
imj = i - j
!xx   10 SUM = SUM + Z(I)*Z(IMJ)                                           
sum = sum + z(i)*z(imj)
10 continue
!xx   20 C(II) = SUM                                                       
c(ii) = sum
20 continue
!                                                                       
!       COVARIANCE MATRIX R COMPUTATION                                 
!                                                                       
!xx      DO 70  II=1,K1                                                    
do 71  ii=1,k1
kmi = k1 - ii + 1
nmi = n1 - ii
do 70  jj=ii,k1
jmi = jj - ii
sum = c(jmi+1)
if( kmi .gt. l )  go to 40
do 30  kk=kmi,l
ki = kk - jmi
!xx   30 SUM = SUM + Z(KK)*Z(KI)                                           
sum = sum + z(kk)*z(ki)
30 continue
40 continue
if( nmlp1 .gt. nmi )  go to 60
do 50  kk=nmlp1,nmi
ki = kk - jmi
!xx   50 SUM = SUM + Z(KK)*Z(KI)                                           
sum = sum + z(kk)*z(ki)
50 continue
60 r(ii,jj) = sum
r(jj,ii) = sum
70 continue
71 continue
!                                                                       
!       DAVIDON'S MINIMIZATION PROCEDURE                                
!                                                                       
f0 = 1.0d60
do  80   i=1,5
!                                                                       
!c      CALL  DAVIDN( FUNCT,HESIAN,A,K,IHES )                             
!cx      CALL  DAVIDN( FUNCT,HESIAN,Z,N,A,K,R,IHES,ISW,IPR,AIC,SD,IFG,LU )
!x      CALL  DAVIDN( FUNCT,HESIAN,Z,N,A,K,R,IHES,ISW,IPR,AIC,SD,
!x     *                     IFG,LU,JER )
!xx      CALL  DAVIDN( FUNCT,HESIAN,Z,N,A,K,R,IHES,ISW,IPR,AIC,SD,JER )
call  davidn( funct,hesian,z,n,a,k,r,ihes,isw,aic,sd,jer )
if( jer .ne. 0 ) return
!                                                                       
if( f0-aic .lt. 0.001 )     go to 90
!xx   80 F0 = AIC
f0 = aic
80 continue
90 continue
!                                                                       
!c      CALL  PRINTA( A,SD,K,TTL,8,TITLE,1,N )                            
sdm = sd
!                                                                       
return
!c    3 FORMAT( //1H ,29(1H-),/,' EXACT LIKELIHOOD MAXIMIZATION',/,1H ,29(
!c     11H-) )                                                            
!xx    3 FORMAT( //1H ,29(1H-),/,' EXACT LIKELIHOOD MAXIMIZATION',10x,
!xx     *'ORDER OF THE AR-MODEL =',i5,/1H ,29(1H-) )                                                            
!xx   65 FORMAT( 1H ,10D13.5 )                                             
end
!c      SUBROUTINE  DAVIDN( FUNCT,HESIAN,X,N,IHES )                       
!xx      SUBROUTINE  DAVIDN( FUNCT,HESIAN,Z,NZ,X,N,R,IHES,ISW,IPR,AIC,SD,
!cx     *                    IFG,LU )
!x     *                    IFG,LU,JER )
!xx     *                    JER )
subroutine  davidn( funct,hesian,z,nz,x,n,r,ihes,isw,aic,sd,jer )
  use timsac_kinds, only: dp
  implicit none
!.....DATE OF THE LATEST REVISION:  JUN. 29, 1979.......................
!                                                                       
!          MINIMIZATION BY DAVIDON-FLETCHER-POWELL PROCEDURE            
!                                                                       
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             FUNCT                                                     
!             HESIAN                                                    
!             LINEAR                                                    
!       ----------------------------------------------------------------
!          INPUTS:                                                      
!             FUNCT:   EXTERNAL FUNCTION SPECIFICATION                  
!             HESIAN:  EXTERNAL FUNCTION SPECIFICATION                  
!             X:       VECTOR OF INITIAL VALUES                         
!             N:       DIMENSION OF THE VECTOR X                        
!             IHES:    =0   INVERSE OF HESSIAN MATRIX IS NOT AVAILABLE  
!                      =1   INVERSE OF HESSIAN MATRIX IS AVAILABLE      
!                                                                       
!          OUTPUT:                                                      
!             X:       VECTOR OF MINIMIZING SOLUTION                    
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  X(50) , DX(50) , G(50) , G0(50) , Y(50)                
!c      DIMENSION  H(50,50) , WRK(50) , S(50)                             
!c      DIMENSION  R(31,31)                                               
!xx      DIMENSION  Z(NZ)
!xx      DIMENSION  X(N) , DX(N) , G(N) , G0(N) , Y(N)                
!xx      DIMENSION  H(N,N) , WRK(N) , S(N)                             
!xx      DIMENSION  R(N+1,N+1)
integer nz, n, ihes, isw, jer
real(dp) z(nz), x(n), r(n+1,n+1), aic, sd
! local
integer i, ic, icc, ig, j
real(dp) dx(n), g(n), g0(n), y(n), h(n,n), wrk(n), s(n),&
&tau1, tau2, eps1, eps2, ramda, const1, sum,&
&s1, s2, stem, ss, ds2, gtem, ed, f, xm, xmb
!
!c      COMMON     / CCC /  ISW, IPR                                      
!c      COMMON     / DDD /  R , F , AIC , SD                              
data  tau1 , tau2  /  1.0d-5 , 1.0d-5  /
data  eps1 , eps2  / 1.0d-5 , 1.0d-5  /
external  funct
!
ramda = 0.5d0
const1 = 1.0d-70
!                                                                       
!          INITIAL ESTIMATE OF INVERSE OF HESSIAN                       
!
h(1:n,1:n) = 0.0d00
s(1:n) = 0.0d00
dx(1:n) = 0.0d00
do  20   i=1,n
!xx      DO  10   J=1,N                                                    
!xx   10 H(I,J) = 0.0D00
!xx      S(I) = 0.0D00                                                     
!xx      DX(I) = 0.0D00                                                    
!xx   20 H(I,I) = 1.0D00                                                   
h(i,i) = 1.0d00
20 continue
isw = 0
!                                                                       
!c      CALL  FUNCT( N,X,XM,G,IG )                                        
!x      CALL  FUNCT( Z,NZ,N,X,R,ISW,XM,G,AIC,SD,F,IG )
call  funct( z,nz,n,x,r,isw,xm,g,aic,sd,f,ig,jer )
if( jer .ne. 0 ) return
!                                                                       
!c      IF( IPR .GE. 2 )   WRITE( 6,340 )   XM, SD, AIC                   
!x      IF( (IPR.GE.2) .AND. (IFG.NE.0) )  WRITE( LU,340 )  XM, SD, AIC
!                                                                       
!          INVERSE OF HESSIAN COMPUTATION (IF AVAILABLE)                
!                                                                       
!c      IF( IHES .EQ. 1 )   CALL  HESIAN( X,N,H )                         
if( ihes .eq. 1 )   call  hesian( x,n,nz,r,sd,h )
!                                                                       
icc = 0
!      ITERATION                                                        
11110 continue
icc = icc + 1
do  11111   ic=1,n
if( ic .eq. 1 .and. icc .eq. 1 )     go to 120
!                                                                       
do  40   i=1,n
!xx   40 Y(I) = G(I) - G0(I)                                               
y(i) = g(i) - g0(i)
40 continue
do  60   i=1,n
sum = 0.0d00
do  50   j=1,n
!xx   50 SUM = SUM + Y(J) * H(I,J)                                         
sum = sum + y(j) * h(i,j)
50 continue
!xx   60 WRK(I) = SUM                                                      
wrk(i) = sum
60 continue
s1 = 0.0d00
s2 = 0.0d00
do  70   i=1,n
s1 = s1 + wrk(i) * y(i)
!xx   70 S2 = S2 + DX(I) * Y(I)                                            
s2 = s2 + dx(i) * y(i)
70 continue
if( s1.le.const1 .or. s2.le.const1 )  go to 900
if( s1 .le. s2 )     go to 100
!                                                                       
!          UPDATE THE INVERSE OF HESSIAN MATRIX                         
!                                                                       
!               ---  DAVIDON-FLETCHER-POWELL TYPE CORRECTION  ---       
!                                                                       
!xx      DO  90   I=1,N                                                    
do  91   i=1,n
do  90   j=i,n
h(i,j) = h(i,j) + dx(i)*dx(j)/s2 - wrk(i)*wrk(j)/s1
!xx   90 H(J,I) = H(I,J)                                                   
h(j,i) = h(i,j)
90 continue
91 continue
go to  120
!                                                                       
!               ---  FLETCHER TYPE CORRECTION  ---                      
!                                                                       
100 continue
stem = s1 / s2 + 1.0d00
!xx      DO  110   I=1,N                                                   
do  111   i=1,n
do  110   j=i,n
h(i,j) = h(i,j)- (dx(i)*wrk(j)+wrk(i)*dx(j)-dx(i)*dx(j)*stem)/s2
!xx  110 H(J,I) = H(I,J)                                                   
h(j,i) = h(i,j)
110 continue
111 continue
!                                                                       
!                                                                       
!                                                                       
120 continue
ss = 0.0d00
do  150   i=1,n
sum = 0.0d00
do  140   j=1,n
!xx  140 SUM = SUM + H(I,J)*G(J)                                           
sum = sum + h(i,j)*g(j)
140 continue
ss = ss + sum * sum
!xx  150 S(I) = -SUM                                                       
s(i) = -sum
150 continue
!                                                                       
!                                                                       
s1 = 0.0d00
s2 = 0.0d00
do  170   i=1,n
s1 = s1 + s(i)*g(i)
!xx  170 S2 = S2 + G(I)*G(I)
s2 = s2 + g(i)*g(i)
170 continue
ds2 = dsqrt(s2)
gtem = dabs(s1) / ds2
if( gtem .le. tau1  .and.  ds2 .le. tau2 )     go to  900
if( s1 .lt. 0.0d00 )     go to  200
h(1:n,1:n) = 0.0d00
do  190   i=1,n
!xx      DO  180   J=1,N                                                   
!xx  180 H(I,J) = 0.0D00                                                   
h(i,i) = 1.0d00
!xx  190 S(I) = -S(I)                                                      
s(i) = -s(i)
190 continue
200 continue
!                                                                       
ed = xm
!                                                                       
!          LINEAR  SEARCH                                               
!                                                                       
!c      CALL  LINEAR( FUNCT,X,S,RAMDA,ED,N,IG )                           
!xx      CALL  LINEAR( FUNCT,Z,NZ,X,S,RAMDA,ED,N,R,AIC,SD,F,ISW,IPR,IG,
!cx     *              IFG,LU )
!x     *              IFG,LU,JER )
!xx     *              JER )
call  linear( funct,z,nz,x,s,ramda,ed,n,r,aic,sd,f,isw,ig,jer )
if( jer .ne. 0 ) return
!                                                                       
!c      IF( IPR .GE. 2 )   WRITE( 6,330 )   RAMDA, F, SD, AIC             
!x      IF( (IPR.GE.2) .AND. (IFG.NE.0) )
!x     *                     WRITE( LU,330 )  RAMDA, F, SD, AIC
!                                                                       
s1 = 0.0d00
do  210   i=1,n
dx(i) = s(i) * ramda
s1 = s1 + dx(i) * dx(i)
g0(i) = g(i)
!xx  210 X(I) = X(I) + DX(I)                                               
x(i) = x(i) + dx(i)
210 continue
xmb = xm
isw = 0
!                                                                       
!c      CALL  FUNCT( N,X,XM,G,IG )                                        
!x      CALL  FUNCT( Z,NZ,N,X,R,ISW,XM,G,AIC,SD,F,IG )
call  funct( z,nz,n,x,r,isw,xm,g,aic,sd,f,ig,jer )
if( jer .ne. 0 ) return
!                                                                       
s2 = 0.d0
do  220     i=1,n
!xx  220 S2 = S2 + G(I)*G(I)
s2 = s2 + g(i)*g(i)
220 continue
if( dsqrt(s2) .gt. tau2 )   go to  11111
if( xmb/xm-1.d0 .lt. eps1  .and.  dsqrt(s1) .lt. eps2 )  go to 900
11111 continue
if( icc .ge. 5 )     go to 900
go to 11110
900 continue
!x      IF( IPR .LE. 0 )   RETURN                                         
!c      WRITE( 6,600 )                                                    
!c      WRITE( 6,610 )     (X(I),I=1,N)                                   
!c      WRITE( 6,601 )                                                    
!c      WRITE( 6,610 )     (G(I),I=1,N)                                   
!x      IF( IFG .EQ. 0 )   RETURN
!x      WRITE( LU,600 )                                                    
!x      WRITE( LU,610 )     (X(I),I=1,N)                                   
!x      WRITE( LU,601 )                                                    
!x      WRITE( LU,610 )     (G(I),I=1,N)                                   
return
!xx  330 FORMAT( 1H ,'LAMBDA =',D15.7,3X,'(-1)LOG LIKELIHOOD =',D23.15,3X, 
!xx     * 'SD =',D22.15,5X,'AIC =',D23.15 )                                
!xx  340 FORMAT( 1H ,26X,'(-1)LOG-LIKELIHOOD =',D23.15,3X,'SD =',D22.15,5X,
!xx     *  'AIC =',D23.15 )                                                
!xx  600 FORMAT( 1H ,'-----  X  -----' )                                   
!xx  601 FORMAT( 1H ,'***  GRADIENT  ***' )                                
!xx  610 FORMAT( 1H ,10D13.5 )                                             
end
!c      SUBROUTINE  FUNCT( M,A,F,G,IFG )                                  
!x      SUBROUTINE  FUNCT( Z,N,M,A,R,ISW,F,G,AIC,SD,FF,IFG )
subroutine  funct( z,n,m,a,r,isw,f,g,aic,sd,ff,ifg,jer )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE COMPUTES THE EXACT LIKELIHOOD AND ITS GRADIENT OF 
!     THE M-TH ORDER AR-MODEL.                                          
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             ARCOEF                                                    
!             PARCOR                                                    
!             SUBDET                                                    
!       ----------------------------------------------------------------
!                                                                       
!          INPUTS:                                                      
!             M:    ORDER OF THE AR MODEL                               
!             A:    VECTOR OF AR-COEFFICIENTS                           
!                                                                       
!          OUTPUTS:                                                     
!             F:    LIKELIHOOD OF THE AR-MODEL                          
!             G:    GRADIENT OF THE LIKELIHOOD FUNCTION                 
!             IFG:  =0    ;IF MODEL IS STATIONALY                       
!                   =1    ;IF MODEL IS NON-STATIONALY                   
!                                                                       
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!C      REAL * 4   Z                                                      
!c      DIMENSION  Z(10000)
!C      COMMON     / AAA /  N , Z(1)                                      
!c      COMMON     / AAA /  N
!c      COMMON     / BBB /  Z
!c      COMMON     / CCC /  ISW, IDUMMY                                          
!c      COMMON     / DDD /  R , FF , AIC , SD                             
!c      DIMENSION   R(31,31) , T(31,31) , U(31,31) , S(31,31)              
!c      DIMENSION   A(20) , B(30) , G(30)                                  
!xx      DIMENSION   Z(N)
!xx      DIMENSION   R(M+1,M+1) , T(M+1,M+1) , U(M+1,M+1) , S(M+1,M+1)              
!xx      DIMENSION   A(M) , B(M) , G(M)                                  
integer n, m, isw, ifg, jer
real(dp) z(n), a(m), r(m+1,m+1), f, g(m), aic, sd, ff
! local
integer i, ii, im1, j, j1, jc, je, je1, jj, jm1, k, kk, km1,&
&m1, mj, mp1
real(dp) t(m+1,m+1), u(m+1,m+1), s(m+1,m+1), b(m),&
&dn, dn1, sum, dsd, usum, dett, udet
!c      MJ = 31                                                           
mj = m+1
!                                                                       
280 ifg = 0
!                                                                       
dn = n
dn1 = 1.0d00 / dn
mp1 = m + 1
m1 = mp1 / 2
!                                                                       
!   **  COMPUTATION OF F  **                                            
!                                                                       
!                                                                       
!          INVERSE OF COVARIANCE MATRIX COMPUTATION                     
!                                                                       
!                                                                       
t(1,1) = 1.0d00 - a(m) * a(m)
if( m .eq. 1 )     go to 45
do  10     i=2,m
im1 = i - 1
ii = m - im1
!xx   10 T(I,1) = - ( A(IM1)+A(M)*A(II) )                                  
t(i,1) = - ( a(im1)+a(m)*a(ii) )
10 continue
if( m .eq. 2 )     go to 25
!                                                                       
!xx      DO  20     J=2,M1                                                 
do  21     j=2,m1
jm1 = j - 1
jj = m - jm1
do  20     i=j,jj
im1 = i - 1
ii = m - im1
!xx   20 T(I,J) = T(IM1,JM1) + A(IM1)*A(JM1) - A(II)*A(JJ)                 
t(i,j) = t(im1,jm1) + a(im1)*a(jm1) - a(ii)*a(jj)
20 continue
21 continue
!                                                                       
!xx   25 DO  30     J=1,M1                                                 
25 do  31     j=1,m1
je = mp1 - j
je1 = je - 1
do  30     i=j,je1
ii = mp1 - i
!xx   30 T(JE,II) = T(I,J)                                                 
t(je,ii) = t(i,j)
30 continue
31 continue
!xx      DO  40     I=2,M                                                  
do  41     i=2,m
im1 = i - 1
do  40     j=1,im1
!xx   40 T(J,I) = T(I,J)                                                   
t(j,i) = t(i,j)
40 continue
41 continue
!                                                                       
45 continue
!                                                                       
!         INNOVATION VARIANCE COMPUTATION                               
!                                                                       
sd = r(1,1)
do  60     j=1,m
j1 = j + 1
sum = 0.0d00
do  50     k=1,m
!xx   50 SUM = SUM + A(K)*R(J1,K+1)                                        
sum = sum + a(k)*r(j1,k+1)
50 continue
!xx   60 SD = SD + SUM*A(J)                                                
sd = sd + sum*a(j)
60 continue
sum = 0.0d00
do  70     k=1,m
!xx   70 SUM = SUM + A(K)*R(1,K+1)                                         
sum = sum + a(k)*r(1,k+1)
70 continue
sd = sd - 2.0d00*sum
!                                                                       
do  80     j=1,m
sum = 0.0d00
do  85     k=1,m
!xx   85 SUM = SUM + T(J,K)*Z(K)                                           
sum = sum + t(j,k)*z(k)
85 continue
!xx   80 SD = SD + SUM*Z(J)                                                
sd = sd + sum*z(j)
80 continue
sd = sd * dn1
!                                                    DS     ) 033,6 (ETI
!xx      DO  90     I=1,M                                                  
do  91     i=1,m
do  90     j=1,m
!xx   90 U(I,J) = T(I,J)                                                   
u(i,j) = t(i,j)
90 continue
91 continue
call  subdet( u,dett,m,mj )
!                                                                       
!         CHECK THE STATIONARITY                                        
!                                                                       
!                                                                       
call  parcor( a,m,b )
!                                                                       
do  210     i=1,m
if( dabs(b(i)) .lt. 1.d0 )     go to  210
if( isw .eq. 0 )   b(i) = 0.999d0*dabs(b(i))/b(i)
ifg = 1
210 continue
if( ifg .eq. 0 )     go to  220
if( isw .eq. 1 )     return
!                                                                       
call  arcoef( b,m,a )
go to  280
220 continue
if( sd .gt. 0.0d00 ) go to 221
jer = 11111
return
221 continue
!                                                                       
!         LIKELIHOOD                                                    
!                                                                       
f = 0.5d00 * ( dn*dlog(sd) - dlog(dett) )
ff = f
aic = 2.d0*f + 2*(m+1)
!                                                                       
if( isw .eq. 1 )     return
!                                                                       
!       COMPUTATION OF GRADIENT OF F                                    
!                                                                       
do  1000     jc=1,m
!                                                                       
s(1,1) = 0.0d00
if( jc .eq. m )     s(1,1) = -2.0d00*a(m)
if( m .eq. 1 )     go to 145
!                                                                       
do  100     i=2,m
ii = mp1 - i
im1 = i - 1
s(i,1) = 0.0d00
if( ii .eq. jc )     s(i,1) = -a(m)
if( m .eq. jc )     s(i,1) = -a(ii)
if( im1 .eq. jc )     s(i,1) = s(i,1) - 1.0d00
100 continue
if( m .eq. 2 )     go to 120
!                                                                       
!xx      DO  110     K=2,M1                                                
do  111     k=2,m1
km1 = k - 1
kk = mp1 - k
do  110     i=k,kk
im1 = i - 1
ii = mp1 - i
s(i,k) = s(im1,km1)
if( km1 .eq. jc )     s(i,k) = s(i,k) + a(im1)
if( im1 .eq. jc )     s(i,k) = s(i,k) + a(km1)
if( ii .eq. jc )     s(i,k) = s(i,k) - a(kk)
if( kk .eq. jc )     s(i,k) = s(i,k) - a(ii)
110 continue
111 continue
!                                                                       
!xx  120 DO  130     J=1,M1                                                
120 do  131     j=1,m1
je=  mp1 - j
je1 = je - 1
do  130     i=j,je1
ii = mp1 - i
!xx  130 S(JE,II) = S(I,J)                                                 
s(je,ii) = s(i,j)
130 continue
131 continue
!xx      DO  140     I=2,M                                                 
do  141     i=2,m
im1 = i - 1
do  140     j=1,im1
!xx  140 S(J,I) = S(I,J)                                                   
s(j,i) = s(i,j)
140 continue
141 continue
145 continue
!                                                                       
!                                                                       
j1 = jc + 1
dsd = r(1,j1)
do  150     k=1,m
!xx  150 DSD = DSD - A(K) * R(J1,K+1)                                      
dsd = dsd - a(k) * r(j1,k+1)
150 continue
!                                                                       
dsd = -2.0d00*dsd
do  170     i=1,m
sum = 0.0d00
do  160     j=1,m
!xx  160 SUM = SUM + S(I,J)*Z(J)
sum = sum + s(i,j)*z(j)
160 continue
!xx  170 DSD = DSD + SUM*Z(I)
dsd = dsd + sum*z(i)
170 continue
dsd = dsd * dn1
usum = 0.0d00
do  200     k=1,m
!xx      DO  180     I=1,M
do  181     i=1,m
do  180     j=1,m
!xx  180 U(I,J) = T(I,J)
u(i,j) = t(i,j)
180 continue
181 continue
do  190     j=1,m
!xx  190 U(K,J) = S(K,J)                                                   
u(k,j) = s(k,j)
190 continue
!                                                                       
call  subdet( u,udet,m,mj )
!                                                                       
!xx  200 USUM = USUM + UDET
usum = usum + udet
200 continue
!                                                                       
!         GRADIENT OF LIKELIHOOD FUNCTION                               
!                                                                       
g(jc) = 0.5d00*( -dn*dsd/sd + usum/dett )
1000 continue
do  300     i=1,m
!xx  300 G(I) = -G(I)
g(i) = -g(i)
300 continue
return
end
!c      SUBROUTINE  HESIAN( X,K,H )                                       
subroutine  hesian( x,k,n,r,sd,h )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE RETURNS THE INVERSE OF AN APPROXIMATION TO THE HES
!     OF LOG-LIKELIHOOD FUNCTION OF THE AUTOREGRESSIVE MODEL OF ORDER K.
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
!             INVDET                                                    
!       ----------------------------------------------------------------
!                                                                       
!       INPUTS:                                                         
!          X:     VECTOR OF AR-COEFFICIENTS                             
!          K:     ORDER OF THE AR-MODEL                                 
!          N:     DATA LENGTH                                           
!          Z:     VECTOR OF ORIGINAL DATA                               
!          R:     SQUARE ROOT OF THE COVARIANCE MATRIX, OUTPUT OF SUBROU
!                 ARMLE                                                 
!          SD:    ESTIMATE OF INNOVATION VARIANCE                       
!                                                                       
!       OUTPUT:                                                         
!          H:     INVERSE OF THE HESSIAN                                
!xx      IMPLICIT  REAL*8( A-H,O-Z )                                       
!C      REAL*4  Z                                                         
!c      DIMENSION  Z(10000)
!C      COMMON     / AAA /  N , Z(1)                                      
!c      COMMON     / AAA /  N
!c      COMMON     / BBB /  Z
!c      COMMON     / DDD /  R , F , AIC , SD                              
!c      DIMENSION  R(31,31)                                               
!x      DIMENSION  X(1)                                                   
!c      DIMENSION  H(50,50) , S(50)                                       
!xx      DIMENSION  R(K+1,K+1)                                               
!xx      DIMENSION  X(K)                                                   
!xx      DIMENSION  H(K,K) , S(K)
integer k, n
real(dp) x(k), r(k+1,k+1), sd, h(k,k)
! local
integer i, ii, j
real(dp) s(k), sum, hdet
do 10  i=1,k
sum = r(1,i+1)
do 20  ii=1,k
!xx   20 SUM = SUM - X(II)*R(II+1,I+1)                                     
sum = sum - x(ii)*r(ii+1,i+1)
20 continue
!xx   10 S(I) = SUM / SD
s(i) = sum / sd
10 continue
!xx      DO 30  I=1,K                                                      
do 31  i=1,k
do 30  j=1,k
!xx   30 H(I,J) = 0.5D0 * (R(I+1,J+1)/SD - S(I)*S(J)/N)                    
h(i,j) = 0.5d0 * (r(i+1,j+1)/sd - s(i)*s(j)/n)
30 continue
31 continue
!c      CALL  INVDET( H,HDET,K,50 )                                       
call  invdet( h,hdet,k,k )
return
end
!c      SUBROUTINE  LINEAR( FUNCT,X,H,RAM,EE,K,IG )                       
!xx      SUBROUTINE  LINEAR( FUNCT,Z,N,X,H,RAM,EE,K,R,AIC,SD,F,ISW,IPR,IG,
!cx     *                    JFG,LU )
!x     *                    JFG,LU,JER )
!xx     *                    
subroutine  linear( funct,z,n,x,h,ram,ee,k,r,aic,sd,f,isw,ig,jer )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE PERFORMS THE LINEAR SEARCH ALONG THE DIRECTION SPE
!     BY THE VECTOR H                                                   
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
!             FUNCT                                                     
!       ----------------------------------------------------------------
!                                                                       
!     INPUTS:                                                           
!        FUNCT:   EXTERNAL FUNCTION SPECIFICATION                       
!        X:       VECTOR OF POSITION                                    
!        H:       SEARCH DIRECTION                                      
!        K:       DIMENSION OF VECTOR X                                 
!                                                                       
!     OUTPUTS:                                                          
!        RAM:     OPTIMAL STEP WIDTH                                    
!        E2:      MINIMUM FUNCTION VALUE                                
!        IG:      ERROR CODE                                            
!                                                                       
!xx      IMPLICIT  REAL  *8 ( A-H,O-Z )                                    
!xx      INTEGER  RETURN,SUB                                               
!c      DIMENSION  X(1) , H(1) , X1(50)                                   
!c      DIMENSION  G(50)                                                  
!c      COMMON     / CCC /  ISW , IPR                                     
!c      DIMENSION  R(31,31)
!xx      DIMENSION  X(K) , H(K) , X1(K)                                   
!xx      DIMENSION  G(K)                                                  
!xx      DIMENSION  Z(N)
!xx      DIMENSION  R(K+1,K+1)

integer n, k, isw, ig, jer
real(dp) z(n), x(k), h(k), ram, ee, r(k+1,k+1), aic, sd, f
! local
integer i, ifg, return, sub
real(dp) x1(k), g(k), const2, hnorm, e1, e2, e3,&
&ram1, ram2, ram3, a1, a2, a3, b1, b2
!                                                                       
isw = 1
if( ram .le. 1.0d-30 )  ram = 0.01d0
const2 = 1.0d-60
hnorm = 0.d0
do 10  i=1,k
!xx   10 HNORM = HNORM + H(I)**2                                           
hnorm = hnorm + h(i)**2
10 continue
hnorm = dsqrt( hnorm )
!                                                                       
ram2 = ram
e1 =ee
ram1 = 0.d0
!                                                                       
do 20  i=1,k
!xx   20 X1(I) = X(I) + RAM2*H(I)                                          
x1(i) = x(i) + ram2*h(i)
20 continue
!c      CALL  FUNCT( K,X1,E2,G,IG )                                       
!c      IF(IPR.GE.7)  WRITE(6,2)  RAM2,E2                                 
!cx      CALL  FUNCT( Z,N,K,X1,R,ISW,E2,G,AIC,SD,F,IG )
call  funct( z,n,k,x1,r,isw,e2,g,aic,sd,f,ig,jer )
if( jer .ne. 0 ) return
!x      IF( (IPR.GE.7) .AND. (JFG.NE.0) )  WRITE(LU,2)  RAM2,E2
!                                                                       
if( ig .eq. 1 )  go to  50
if( e2 .gt. e1 )  go to 50
30 ram3 = ram2*2.d0
do 40  i=1,k
!xx   40 X1(I) = X(I) + RAM3*H(I)                                          
x1(i) = x(i) + ram3*h(i)
40 continue
!c      CALL  FUNCT( K,X1,E3,G,IG )                                       
!cx      CALL  FUNCT( Z,N,K,X1,R,ISW,E3,G,AIC,SD,F,IG )
call  funct( z,n,k,x1,r,isw,e3,g,aic,sd,f,ig,jer )
if( jer .ne. 0 ) return
if( ig.eq.1 )  go to  500
!c      IF( IPR.GE.7 )  WRITE(6,3)  RAM3,E3                               
!x      IF( (IPR.GE.7) .AND. (JFG.NE.0 ) )  WRITE(LU,3)  RAM3,E3
if( e3 .gt. e2 )  go to 70
ram1 = ram2
ram2 = ram3
e1 = e2
e2 = e3
go to 30
!                                                                       
50 ram3 = ram2
e3 = e2
ram2 = ram3*0.1d0
if( ram2*hnorm .lt. const2 )  go to  400
do 60  i=1,k
!xx   60 X1(I) = X(I) + RAM2*H(I)
x1(i) = x(i) + ram2*h(i)
60 continue
!c      CALL  FUNCT( K,X1,E2,G,IG )                                       
!c      IF(IPR.GE.7)  WRITE(6,4)  RAM2,E2                                 
!cx      CALL  FUNCT( Z,N,K,X1,R,ISW,E2,G,AIC,SD,F,IG )
call  funct( z,n,k,x1,r,isw,e2,g,aic,sd,f,ig,jer )
if( jer .ne. 0 ) return
!x      IF( (IPR.GE.7) .AND. (JFG.NE.0) )  WRITE(LU,4)  RAM2,E2
if( e2.gt.e1 )  go to 50
!                                                                       
!c   70 ASSIGN 80 TO RETURN                                               
70 return = 80
go to 200
!                                                                       
80 do 90  i=1,k
!xx   90 X1(I) = X(I) + RAM*H(I)                                           
x1(i) = x(i) + ram*h(i)
90 continue
!c      CALL  FUNCT( K,X1,EE,G,IG )                                       
!c      IF(IPR.GE.7)  WRITE(6,5)  RAM,EE                                  
!cx      CALL  FUNCT( Z,N,K,X1,R,ISW,EE,G,AIC,SD,F,IG )
call  funct( z,n,k,x1,r,isw,ee,g,aic,sd,f,ig,jer )
if( jer .ne. 0 ) return
!x      IF( (IPR.GE.7) .AND. (JFG.NE.0) )  WRITE(LU,5)  RAM,EE
!                                                                       
ifg = 0
!c      ASSIGN  300 TO  SUB                                               
!c      ASSIGN 200 TO SUB                                                 
!c   95 ASSIGN 130 TO RETURN                                              
sub = 300
sub = 200
95 return = 130
if( ram .gt. ram2 )  go to 110
if( ee .ge. e2 )  go to 100
ram3 = ram2
ram2 = ram
e3 =e2
e2 =ee
!c      GO TO  SUB,( 200,300 )                                            
if( sub .eq. 200 ) go to  200
if( sub .eq. 300 ) go to  300
!                                                                       
100 ram1 = ram
e1 = ee
!c      GO TO  SUB,( 200,300 )                                            
if( sub .eq. 200 ) go to  200
if( sub .eq. 300 ) go to  300
!                                                                       
110 if( ee .le. e2 )  go to 120
ram3 = ram
e3 = ee
!c      GO TO  SUB,( 200,300 )                                            
if( sub .eq. 200 ) go to  200
if( sub .eq. 300 ) go to  300
!                                                                       
120 ram1 = ram2
ram2 = ram
e1 = e2
e2 = ee
!c      GO TO  SUB,( 200,300 )                                            
if( sub .eq. 200 ) go to  200
if( sub .eq. 300 ) go to  300
!                                                                       
130 do 140  i=1,k
!xx  140 X1(I) = X(I) + RAM*H(I)                                           
x1(i) = x(i) + ram*h(i)
140 continue
!c      CALL  FUNCT( K,X1,EE,G,IG )                                       
!c      IF( IPR.GE.7 )  WRITE(6,6)  RAM,EE                                
!cx      CALL  FUNCT( Z,N,K,X1,R,ISW,EE,G,AIC,SD,F,IG )
call  funct( z,n,k,x1,r,isw,ee,g,aic,sd,f,ig,jer )
if ( jer .ne. 0 ) return
!x      IF( (IPR.GE.7)  .AND. (JFG.NE.0) )  WRITE(LU,6)  RAM,EE
!c      ASSIGN 200 TO SUB                                                 
sub = 200
ifg = ifg+1
ifg = 0
if( ifg .eq. 1 )  go to 95
!                                                                       
if( e2 .lt. ee )  ram = ram2
return
!                                                                       
!      -------  INTERNAL SUBROUTINE SUB1  -------                       
200 a1 = (ram3-ram2)*e1
a2 = (ram1-ram3)*e2
a3 = (ram2-ram1)*e3
b2 = (a1+a2+a3)*2.d0
b1 = a1*(ram3+ram2) + a2*(ram1+ram3) + a3*(ram2+ram1)
if( b2 .eq. 0.d0 )  go to 210
ram = b1 /b2
!c      GO TO RETURN ,( 80,130 )                                          
if( return .eq. 80 ) go to 80
if( return .eq. 130 ) go to 130
!                                                                       
210 ig = 1
ram = ram2
return
!                                                                       
!      -------  INTERNAL SUBROUTINE SUB2  -------                       
!                                                                       
300 if( ram3-ram2 .gt. ram2-ram1 )  go to 310
ram = (ram1+ram2)*0.5d0
!c      GO TO RETURN ,( 80,130 )                                          
if( return .eq. 80 ) go to 80
if( return .eq. 130 ) go to 130
!                                                                       
310 ram = (ram2+ram3)*0.5d0
!c      GO TO RETURN ,( 80,130 )                                          
if( return .eq. 80 ) go to 80
if( return .eq. 130 ) go to 130
! ------------------------------------------------------------          
!                                                                       
400 ram = 0.d0
return
! ------------------------------------------------------------          
!                                                                       
500 ram = (ram2+ram3)*0.5d0
510 do 520  i=1,k
!xx  520 X1(I) = X(I) + RAM*H(I)
x1(i) = x(i) + ram*h(i)
520 continue
!c      CALL  FUNCT( K,X1,E3,G,IG )                                       
!c      IF( IPR.GE.7 )  WRITE(6,7)  RAM,E3                                
!cx      CALL  FUNCT( Z,N,K,X1,R,ISW,E3,G,AIC,SD,F,IG )
call  funct( z,n,k,x1,r,isw,e3,g,aic,sd,f,ig,jer )
if( jer .ne. 0 ) return
!x      IF( (IPR.GE.7) .AND. (JFG.NE.0) )  WRITE(LU,7)  RAM,E3
if( ig.eq.1 )  go to 540
if( e3.gt.e2 )  go to 530
ram1 = ram2
ram2 = ram
e1 = e2
e2 = e3
go to 500
!                                                                       
530 ram3 = ram
go to 70
!                                                                       
540 ram = (ram2+ram)*0.5d0
go to 510
!                                                                       
! ------------------------------------------------------------          
!xx    1 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E1 =',D25.17 )                
!xx    2 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E2 =',D25.17 )                
!xx    3 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E3 =',D25.17 )                
!xx    4 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E4 =',D25.17 )                
!xx    5 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E5 =',D25.17 )                
!xx    6 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E6 =',D25.17 )                
!xx    7 FORMAT( 1H ,'LAMBDA =',D18.10, 10X,'E7 =',D25.17 )                
end
