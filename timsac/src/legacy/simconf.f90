! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine simconf(d,k,h,l,r,b0,bx0,s0,q0,bc,bd,g,avy,si,s2)
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM SIMCON                                                    
!     PROGRAM 74.3.2.  OPTIMAL CONTROLLER DESIGN AND SIMULATION         
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO                                                         
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977                    
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN                       
!        "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(2
!        BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
!        NO.6 MARCH 1976, THE INSTITUTE OF STATISTICAL MATHEMATICS      
!-----------------------------------------------------------------------
!     THIS PROGRAM PRODUCES OPTIMAL CONTROLLER GAIN AND SIMULATES THE   
!     CONTROLLED PROCESS. THE BASIC STATE SPACE MODEL IS OBTAINED       
!     FROM THE AUTOREGRESSIVE MOVING AVERAGE MODEL OF A VECTOR          
!     PROCESS Y(I); Y(I)+B(1)Y(I-1)+..+B(K)Y(I-K)=X(I)+A(1)X(I-1)+..    
!     ..+A(K-1)X(I-K+1).                                                
!                                                                       
!     THE FOLLOWING INPUTS ARE REQUIRED:                                
!      THE OUTPUTS OF PROGRAM MARKOV:                                   
!       (D,K): INTEGERS                                                 
!             D, DIMENSION OF Y(I) (LESS THAN OR EQUAL TO MJD)          
!             K, ORDER OF THE PROCESS (LESS THAN OR EQUAL TO MJK)       
!       (B(I),I=1,K): MATRICES OF AUTOREGRESSIVE COEFFICIENTS           
!       (W(I),I=1,K-1): IMPULSE RESPONCE MATRICES OF THE VECTOR VARIABLE
!                       Y(I) TO THE INNOVATION INPUT                    
!      S: D*D COVARIANCE MATRIX OF INNOVATION X(I)                      
!      (H,L,R): INTEGERS                                                
!              H, SPAN OF CONTROL PERFORMANCE EVALUATION                
!              L, LENGTH OF EXPERIMENTAL OBSERVATION                    
!              R, DIMENSION OF CONTROL INPUT (LESS THAN OR EQUAL TO D)  
!                 THE LAST R COMPONENTS OF Y(I) REPRESENT THE CONTROL IN
!                 AND THE REST REPRESENT THE OUTPUT OF THE SYSTEM.      
!      QX: WEIGHTING MATRIX OF PERFORMANCE, POSITIVE DEFINITE           
!                                                                       
!     BASIC EQUATION OF THE SYSTEM IS AS FOLLOWS;                       
!     V(I+1)=A*V(I)+BD*SI(I)+BC*CI(I+1): STATE VECTOR                   
!     CI(I+1)=G*V(I)                   : CONTROLLER INPUT               
!     Y(I)=C*V(I)                      : OBSERVED OUTPUT                
!     WHERE        +-                              -+                   
!              A = I     0       I       0 ..     0 I ,                 
!                  I     0       0       I ..     0 I                   
!                  I     .       .       . ..     . I                   
!                  I -B(K) -B(K-1) -B(K-2) .. -B(1) I                   
!                  +-                              -+                   
!                                                                       
!                  +-                -+                                 
!              C = I   I 0 0 ... 0    I,                                
!                  +-                -+                                 
!                                                                       
!              +-     -+     +-       -+                                
!              I       I     I  I      I                                
!              I       I     I  W(1)   I                                
!              I BD BC I  =  I  W(2)   I  =  BX,                        
!              I       I     I   .     I                                
!              I       I     I  W(K-1) I                                
!              +-     -+     +-       -+                                
!                                                                       
!              G = INVERSE(-(BC'*Q(H)*BC))*BC'*Q(H)*A,                  
!              Z(I+1) = Q(I)-Q(I)*BC*INVERSE(BC'*Q(I)*BC)*BC'*Q(I) (I=0,
!              Q(I+1) = A'*Z(I+1)*A+Q(0) (I=0,H-1),                     
!                     +-       -+                                       
!              Q(0) = I QX 0..0 I                                       
!                     I  0 0..0 I                                       
!                     I  0 0..0 I                                       
!                     I  . ...0 I                                       
!                     I  0 0..0 I                                       
!                     +-       -+                                       
!              SI(I)= OUTPUT (FIRST D-R COMPONENTS OF Y(I)) INNOVATION. 
!     B(I)(I=1,K) ARE THE AUTOREGRESSIVE COEFFICIENTS OF THE AR-MA      
!     REPRESENTATION OF Y(I).                                           
!     W(J)(J=0,K-1)(W(0)=I) ARE THE IMPULSE RESPONSE MATRICES OF Y(I+J) 
!     THE INNOVATION X(I).                                              
!     BY MODIFYING THE INPUT LOADING PROGRAM, Q(0) MAY BE PUT EQUAL TO A
!     ARBITRARY KD*KD POSITIVE DEFINITE MATRIX.                         
!     THE DIMENSIONS OF THE MATRICES ARE AS FOLLOWES:                   
!     A:  (K*D)*(K*D) MATRIX                                            
!     B:  D*D         MATRIX                                            
!     BD: (K*D)*(D-R) MATRIX                                            
!     BC: (K*D)*R     MATRIX                                            
!     W:  D*D         MATRIX                                            
!     Q:  (K*D)*(K*D) MATRIX                                            
!     QX: D*D         MATRIX                                            
!     Z:  (K*D)*(K*D) MATRIX                                            
!     G:  R*(K*D)     MATRIX                                            
!     S:  D*D         MATRIX                                            
!     V:  (K*D)*1     VECTOR                                            
!     SI: (D-R)*1     VECTOR                                            
!     CI: R*1         VECTOR                                            
!     Y:  D*1         VECTOR                                            
!
!xx      IMPLICIT REAL*8 (A-H,O-Z)                                         
!c      REAL*4 RANDOM                                                     
!xx      REAL*8 RANDM 
!xx      INTEGER D,H,R,H1                                                  
!c      DIMENSION B(10,10,10),BD(100,10),BC(100,10)                       
!c      DIMENSION Q(100,100),Q0(100,100),Z(100,100),G(10,100)             
!c      DIMENSION W1(10)                                                  
!c      DIMENSION S(10,10),WIN(10,10),Y(10,1000),W2(100,100)              
!c      DIMENSION SI(10),CI(10),V(100),AV(100),BCN(100),BDN(100)          
!c      DIMENSION W3(10,10),W4(10,10),S2(10)                              
!c      DIMENSION R1(10,100),R3(10,100)                                   
!c      DIMENSION BX(100,10)                                              
!xx      DIMENSION B0(D,D,K),B(D,D,K),BD(K*D,D-R),BC(K*D,R)
!xx      DIMENSION Q(K*D,K*D),Q0(D,D),Z(K*D,K*D),G(R,K*D)
!xx      DIMENSION W1(D)
!xx      DIMENSION S0(D,D),S(D-R,D-R)
!xx      DIMENSION WIN1(D,D),WIN2(R,D),Y(D,L),W2(K*D,K*D)
!xx      DIMENSION SI(D),CI(R),V(K*D),AV(K*D),AVY(D),BCN(K*D),BDN(K*D)
!xx      DIMENSION W31(R,R),W32(D,D),W41(D,D),W42(R,D),S2(D)
!xx      DIMENSION R1(R,K*D),R3(R,K*D)
!xx      DIMENSION BX0((K-1)*D,D),BX(K*D,D)
integer d, k, h, l, r
real(dp) b0(d,d,k), bx0((k-1)*d,d), s0(d,d), q0(d,d),&
&bc(k*d,r), bd(k*d,d-r), g(r,k*d), avy(d),&
&si(d), s2(d)
! local
integer i, i1, i2, ii1, im1, isw, ix, iy, iyy, iz, j, j1, j2, jrn,&
&jy, jyy, k1, k2, kd, kdr, kix, kiz, kjx, kjy, kjz, kk1,&
&kk2, kk3, kk4, kx, h1
real(dp) randm, b(d,d,k), q(k*d,k*d),z(k*d,k*d), w1(d),&
&s(d-r,d-r), win1(d,d), win2(r,d), y(d,l),&
&w2(k*d,k*d), ci(r), v(k*d), av(k*d), bcn(k*d),&
&bdn(k*d), w31(r,r), w32(d,d), w41(d,d), w42(r,d),&
&r1(r,k*d), r3(r,k*d), bx(k*d,d), cst0, cst1,&
&cstm6, bdmn, al, sum, wdet
!c      EQUIVALENCE(Q0(1,1),BCN(1))                                       
!c      EQUIVALENCE(WIN(1,1),W1(1),SI(1),CI(1))                           
!c      EQUIVALENCE(V(1),S2(1))                                           
!c      EQUIVALENCE(Q(1,1),Z(1,1),Y(1,1))                                 
!c      EQUIVALENCE (R1(1,1),G(1,1))                                      
!c      EQUIVALENCE(BX(1,1),R3(1,1),AV(1))                                
!c      DATA B/1000*0.0D-00/,BC/1000*0.0D-00/,BD/1000*0.0D-00/            
!c      DATA Q0/10000*0.0D-00/,Q/10000*0.0D-00/,Z/10000*0.0D-00/          
!c      DATA S/100*0.0D-00/,SI/10*0.0D-00/,CI/10*0.0D-00/                 
!c      DATA W1/10*0.0D-00/,W2/10000*0.0D-00/                             
!c      DATA R1/1000*0.0D-00/                                             
!c      DATA BX/1000*0.0D-00/,WIN/100*0.0D-00/,W3/100*0.0D-00/            
!     ABSOLUTE DIMENSIONS                                               
!c      MJD=10                                                            
!c      MJK=10                                                            
!c      MJKD=MJK*MJD                                                      
!
!     INPUT / OUTPUT DATA FILE OPEN
!c	CALL SETWND
!c	CALL FLOPN2(NFL)
!c	IF (NFL.EQ.0) GO TO 999
!                                                                       
!     INITIAL CONDITION INPUT                                           
!c      READ(5,800) D,K                                                   
!c      WRITE(6,900)                                                      
kd=k*d
kx=kd-d
!     AR-COEFFICIENT MATRICES B(I)(I=1,K) INPUT                         
isw=1
do  107 ix=1,k
!     COMMON SUBROUTINE CALL                                            
!c      CALL REMATX(WIN,D,D,ISW,MJD,MJD)                                  
do  106 i=1,d
do  105 j=1,d
!c  105 B(I,J,IX)=-WIN(I,J)                                               
!xx  105 B(I,J,IX)=-B0(I,J,IX)
b(i,j,ix)=-b0(i,j,ix)
105 continue
106 continue
107 continue
!                                                                       
cst0=0.0d-00
cst1=1.0d-00
cstm6=-6.0d-00
!
!     initialize data area
!
q(1:kd,1:kd) = cst0
z(1:kd,1:kd) = cst0
si(1:d) = cst0
ci(1:r) = cst0
w1(1:d) = cst0
w2(1:kd,1:kd) = cst0
r1(1:r,1:kd) = cst0
win1(1:d,1:d) = cst0
win2(1:r,1:d) = cst0
w31(1:r,1:r) = cst0
w32(1:d,1:d) = cst0
bx(1:kd,1:d) = cst0
!
!     MATRIX BX SET                                                     
k1=k-1
k2=k1*d
do 102 i=1,d
!xx  102 BX(I,I)=CST1
bx(i,i)=cst1
102 continue
!                                                                       
!     IMPULSE RESPONSE MATRICES W(I),I=1,K-1 INPUT                      
!c      I1=0                                                              
!c      DO  117 IX=1,K1                                                   
!     COMMON SUBROUTINE CALL                                            
!     MATRIX INPUT, ROWWISE                                             
!c      CALL REMATX (WIN,D,D,ISW,MJD,MJD)                                 
!c      I1=I1+D                                                           
!c      DO  116 I=1,D                                                     
!c      I2=I1+I                                                           
!c      DO 115 J=1,D                                                      
!c  115 BX(I2,J)=WIN(I,J)                                                 
!c  116 CONTINUE                                                          
!c  117 CONTINUE                                                          
!xx      DO 116 I=1,(K-1)*D
do 126 i=1,(k-1)*d
do 116 j=1,d
bx(i+d,j)=bx0(i,j)
116 continue
126 continue

!                                                                       
!     INNOVATION COVARIANCE MATRIX INPUT                                
!     COMMON SUBROUTINE CALL                                            
!c      CALL REMATX (S,D,D,ISW,MJD,MJD)                                   
!xx      DO 117 I=1,D-R
do 127 i=1,d-r
do 117 j=1,d-r
s(i,j)=s0(i,j)
117 continue
127 continue
!c      READ (5,800) H,L,R                                                
!c      WRITE(6,901) D,K,H,L,R                                            
!     MATRIX Q(0) SET UP                                                
!     COMMON SUBROUTINE CALL                                            
!     FOR A GENERAL Q(0) REPLACE THE FOLLWING STATEMENT 108-109 BY      
!     CALL REMATX(Q0,KD,KD,ISW,MJKD,MJKD)                               
!c  108 CALL REMATX(W3,D,D,ISW,MJD,MJD)                                   
!                                                                       
!                                                                       
!c      DO  109  I=1,D                                                    
!c      DO  118  J=1,D                                                    
!c  118 Q0(I,J)=W3(I,J)                                                   
!c  109 CONTINUE                                                          
!c      DO 119 I=1,KD
!c      DO 119 J=1,KD
!xx      DO 119 I=1,D
do 120 i=1,d
do 119 j=1,d
!xx  119 Q(I,J)=Q0(I,J)
q(i,j)=q0(i,j)
119 continue
120 continue
!                                                                       
!     MATRIX BC AND BD SET UP                                           
kdr=d-r
do 1120 i=1,kd
do 1119 j=1,kdr
!xx 1119 BD(I,J)=BX(I,J)
bd(i,j)=bx(i,j)
1119 continue
1120 continue
do 1122 i=1,kd
do 1121 j=1,r
j1=kdr+j
!xx 1121 BC(I,J)=BX(I,J1) 
bc(i,j)=bx(i,j1)
1121 continue
1122 continue
!                                                                       
!     INITIAL MATRIX PRINT OUT FOR DEBUGGING                            
!c      DO  707 IX=1,K                                                    
!c      DO  708 I=1,D                                                     
!c      DO  709 J=1,D                                                     
!c  709 WIN(I,J)=-B(I,J,IX)                                               
!c  708 CONTINUE                                                          
!c      WRITE (6,700) IX                                                  
!c      CALL SUBMPR (WIN,D,D,MJD,MJD)                                     
!c  707 CONTINUE                                                          
!c      WRITE(6,702)                                                      
!c      CALL SUBMPR (BC,KD,R,MJKD,MJD)                                    
!c      WRITE(6,703)                                                      
!c      CALL SUBMPR (BD,KD,KDR,MJKD,MJD)                                  
!c      WRITE(6,704)                                                      
!c      CALL SUBMPR(W3,D,D,MJD,MJD)                                       
!c      WRITE(6,706)                                                      
!c      CALL SUBMPR(S,D,D,MJD,MJD)                                        
!     INITIAL PRINT FOR DEBUGGING END                                   
!     ITERATIVE COMPUTATION OF G                                        
h1=h+1
do 200 ix=1,h1
!     R1=BC'*Q R1 STORE                                                 
!c      CALL MULTRX(BC,KD,R,MJKD,MJD,Q,KD,KD,MJKD,MJKD,R1,R,KD,MJD,MJKD,2)
call multrx(bc,kd,r,q,kd,kd,r1,r,kd,2)
!     W3=R1*BC                                                          
!c      CALL MULTRX(R1,R,KD,MJD,MJKD,BC,KD,R,MJKD,MJD,W3,R,R,MJD,MJD,1)   
call multrx(r1,r,kd,bc,kd,r,w31,r,r,1)
!     W3=(INVERSE OF W3)                                                
!c      CALL INVDET(W3,WDET,R,MJD)                                        
call invdets(w31,wdet,r)
!     R3=W3*R1 R3 STORE                                                 
!c      CALL MULTRX(W3,R,R,MJD,MJD,R1,R,KD,MJD,MJKD,R3,R,KD,MJD,MJKD,1)   
call multrx(w31,r,r,r1,r,kd,r3,r,kd,1)
if(ix.eq.h1) go to 210
!     W2=R1'*R3                                                         
!c      CALL MULTRX(R1,R,KD,MJD,MJKD,R3,R,KD,MJD,MJKD,W2,KD,KD,MJKD,MJKD,2
!c     A)                                                                 
call multrx(r1,r,kd,r3,r,kd,w2,kd,kd,2)
!     Z=Q-W2                                                            
!c      CALL SUBTAC(Q,W2,Z,KD,KD,MJKD,MJKD)                               
call subtac(q,w2,z,kd,kd)
!                                                                       
!     W2=A'*Z COMPUTATION:                                              
!                                                                       
!     +-                               -+     +-                        
!     I     0  0  0  0  ..  0 -B(K)'    I     I     Z(1,1) .. Z(1,K)    
!     I     I  0  0  0  ..  0 -B(K-1)'  I     I      .          .       
!     I     0  I  0  0  ..  0 -B(K-2)'  I     I      .          .       
!     I     0  0  I  0  ..  0 -B(K-3)'  I  *  I      .          .       
!     I     0  0  0  I  ..  0 -B(K-4)'  I     I      .          .       
!     I     .  .  .  .  ..  .    .      I     I      .          .       
!     I     0  0  0  0  ..  I -B(1)'    I     I     Z(K,1) .. Z(K,K)    
!     +-                               -+     +-                        
!                                                                       
!         +-                                                    -+      
!         I   -B(K)'Z(K,1)          .... -B(K)'Z(K,K)            I      
!         I   -B(K-1)'Z(K,1)+Z(1,1) .... -B(K-1)'Z(K,K)+Z(1,K)   I      
!      =  I   -B(K-2)'Z(K,1)+Z(2,1) .... -B(K-2)'Z(K,K)+Z(2,K)   I ,    
!         I   -B(K-3)'Z(K,1)+Z(3,1) .... -B(K-3)'Z(K,K)+Z(3,K)   I      
!         I             .           ....           .             I      
!         I             .           ....           .             I      
!         I             .           ....           .             I      
!         I   -B(1)'Z(K,1)+Z(K-1,1) .... -B(1)'Z(K,K)+Z(K-1,K)   I      
!         +-                                                    -+      
!                                                                       
!     WHERE Z(J,L)(J=1,K;L=1,K) ARE THE D*D SUB MATRICES OF THE         
!     BLOCK MATRIX Z.                                                   
!     A'*Z IS STORED IN W2.                                             
iz=k+1
!c      KX=KD-D                                                           
j1=d+1
kiz=0
do  408 iy=1,k
iz=iz-1
do  402 i=1,d
do  401 j=1,d
!c  401 WIN(I,J)=B(I,J,IZ)                                                
!xx  401 WIN1(I,J)=B(I,J,IZ)
win1(i,j)=b(i,j,iz)
401 continue
402 continue
kjz=0
do  407 jy=1,k
do  404 i=1,d
kix=kx+i
do  403 j=1,d
kjx=kjz+j
!c  403 W3(I,J)=Z(KIX,KJX)                                                
!xx  403 W32(I,J)=Z(KIX,KJX)
w32(i,j)=z(kix,kjx)
403 continue
404 continue
!c      CALL MULTRX (WIN,D,D,MJD,MJD,W3,D,D,MJD,MJD,                      
!c     AW4,D,D,MJD,MJD,2)                                                 
call multrx (win1,d,d,w32,d,d,w41,d,d,2)
do  406 i=1,d
kix=kiz+i
do  405 j=1,d
kjx=kjz+j
!c  405 W2(KIX,KJX)=W4(I,J)                                               
!xx  405 W2(KIX,KJX)=W41(I,J)
w2(kix,kjx)=w41(i,j)
405 continue
406 continue
kjz=kjz+d
407 continue
kiz=kiz+d
408 continue
do  411 iy=j1,kd
do  410 jy=1,kd
iyy=iy-d
!xx  409 W2(IY,JY) = W2(IY,JY) + Z(IYY,JY)                                 
w2(iy,jy) = w2(iy,jy) + z(iyy,jy)
410 continue
411 continue
!                                                                       
!     Q=W2*A COMPUTATIONF W2=A'*Z                                       
!                                                                       
!     +-                   -+     +-                                -+  
!     I  W2(1,1)...W2(1,K)  I     I     0     I        0   ..    0   I  
!     I     .          .    I  *  I     0     0        I   ..    0   I  
!     I     .          .    I     I     .     .        .   ..    .   I  
!     I  W2(K,1)...W2(K,K)  I     I  -B(K) -B(K-1) -B(K-2) .. -B(1)  I  
!     +-                   -+     +-                                -+  
!                                                                       
!         +-                                                            
!         I  -W2(1,K)B(K) -W2(1,K)B(K-1)+W2(1,1) ... -W2(1,K)B(1)+W2(1,K
!         I  -W2(2,K)B(K) -W2(2,K)B(K-1)+W2(2,1) ... -W2(2,K)B(1)+W2(2,K
!         I  -W2(3,K)B(K) -W2(3,K)B(K-1)+W2(3,1) ... -W2(3,K)B(1)+W2(3,K
!      =  I  -W2(4,K)B(K) -W2(4,K)B(K-1)+W2(4,1) ... -W2(4,K)B(1)+W2(4,K
!         I         .                  .         ...        .           
!         I         .                  .         ...        .           
!         I         .                  .         ...        .           
!         I  -W2(K,K)B(K) -W2(K,K)B(K-1)+W2(K,1) ... -W2(K,K)B(1)+W2(K,K
!         +-                                                            
!                                                                       
!     WHERE W2(J,L)(J=1,K;L=1,K) ARE THE (J,L) D*D SUB MATRICES         
!     OF THE BLOCK MATRIX W2.                                           
!     W2*A IS STORED IN Q.                                              
iz=k+1
kjz=0
do  419 jy=1,k
iz=iz-1
do  413 i=1,d
do  412 j=1,d
!xx  412 W32(I,J)=B(I,J,IZ)
w32(i,j)=b(i,j,iz)
412 continue
413 continue
kiz=0
do  418 iy=1,k
do  415 j=1,d
kjx=kx+j
do  414 i=1,d
kix=kiz+i
!c  414 WIN(I,J)=W2(KIX,KJX)                                              
!xx  414 WIN1(I,J)=W2(KIX,KJX)
win1(i,j)=w2(kix,kjx)
414 continue
415 continue
!c      CALL MULTRX (WIN,D,D,MJD,MJD,W3,D,D,MJD,MJD,                      
!c     AW4,D,D,MJD,MJD,1)                                                 
call multrx (win1,d,d,w32,d,d,w41,d,d,1)
do  417 j=1,d
kjx=kjz+j
do  416 i=1,d
kix=kiz+i
!c  416 Q(KIX,KJX)=W4(I,J)                                                
!xx  416 Q(KIX,KJX)=W41(I,J)                                                
q(kix,kjx)=w41(i,j)
416 continue
417 continue
kiz=kiz+d
418 continue
kjz=kjz+d
419 continue
do  422 jy=j1,kd
do  421 iy=1,kd
jyy=jy-d
!xx  420 Q(IY,JY)=Q(IY,JY)+W2(IY,JYY)                                      
q(iy,jy)=q(iy,jy)+w2(iy,jyy)
421 continue
422 continue
!c      DO 142 I=1,KD                                                     
!c      DO 142 J=1,KD                                                     
!xx      DO 142 I=1,D
do 143 i=1,d
do 142 j=1,d
!xx  142 Q(I,J)=Q(I,J)+Q0(I,J)
q(i,j)=q(i,j)+q0(i,j)
142 continue
143 continue
200 continue
!                                                                       
!     G=-(INVERSE(BC'*Q*BC))*BC'*Q*A                                    
!     G=R3*A; R3=INVERSE(BC'*Q*BC)*BC'*Q:                               
!                                                                       
!     +-                   -+     +-                               -+   
!     I  W(1) W(2) .. W(K)  I  *  I    0      I       0   ..    0   I   
!     +-                   -+     I    0      0       I   ..    0   I   
!                                 I    0      0       0   ..    0   I   
!                                 I    .      .       .   ..    .   I   
!                                 I -B(K) -B(K-1) -B(K-2) .. -B(1)  I   
!                                 +-                               -+   
!                                                                       
!       +-                                                              
!       I                                                               
!     = I -W(K)B(K) -W(K)B(K-1)+W(1) -W(K)B(K-2)+W(2) ... -W(K)B(1)+W(K-
!       I                                                               
!       +-                                                              
!                                                                       
!     WHERE W(J)(J=1,K) ARE THE J-TH R*D SUB MATRIX OF THE BLOCK MATRIX 
!     R3*A IS STORED IN G.                                              
210 continue
do  455 i=1,r
do  454 j=1,d
kjx=kx+j
!c  454 WIN(I,J)=R3(I,KJX)                                                
!xx  454 WIN2(I,J)=R3(I,KJX)
win2(i,j)=r3(i,kjx)
454 continue
455 continue
iz=k+1
do  461 jy=1,k
iz=iz-1
do  453 i=1,d
do  452 j=1,d
!c  452 W3(I,J)=B(I,J,IZ)                                                 
!xx  452 W32(I,J)=B(I,J,IZ)
w32(i,j)=b(i,j,iz)
452 continue
453 continue
!c      CALL MULTRX (WIN,R,D,MJD,MJD,W3,D,D,MJD,MJD,                      
!c     AW4,R,D,MJD,MJD,1)                                                 
call multrx (win2,r,d,w32,d,d,w42,r,d,1)
kjz=(jy-1)*d
do 457 i=1,r
do 456 j=1,d
kjy=kjz+j
!c  456 G(I,KJY)=W4(I,J)                                                  
!xx  456 G(I,KJY)=W42(I,J)
g(i,kjy)=w42(i,j)
456 continue
457 continue
!c      IF  (KJZ) 461,461,458                                             
if  (kjz .le. 0) go to 461
if  (kjz .gt. 0) go to 458
458 do  460 i=1,r
do  459 j=1,d
kjx=kjz+j-d
kjy=kjz+j
!xx  459 G(I,KJY)=G(I,KJY)+R3(I,KJX)
g(i,kjy)=g(i,kjy)+r3(i,kjx)
459 continue
460 continue
461 continue
do 202 i=1,r
do 201 j=1,kd
!xx  201 G(I,J)=-G(I,J)
g(i,j)=-g(i,j)
201 continue
202 continue
!     G PRINT OUT                                                       
!c      WRITE(6,903) R,KD                                                 
!     COMMON SUBROUTINE CALL                                            
!c      CALL SUBMPR(G,R,KD,MJD,MJKD)                                      
!                                                                       
!                                                                       
!                                                                       
!     ****************                                                  
!     SIMULATION START                                                  
!     ****************                                                  
!     V(I),I=1,KD  CLEAR                                                
!xx      DO  301 I=1,KD                                                    
!xx  301 V(I)=CST0
v(1:kd)=cst0
!     CONSTANT FOR NOISE GENERATION                                     
!     COMMON SUBROUTINE CALL                                            
!c      SI(1)=RANDOM(1)                                                   
si(1)=randm(1,kk1,kk2,kk3,kk4)
!     OBSERVED SYSTEM OUTPUT GENERATION I=1,L                           
!     COVARIANCE MATRIX FACTORIZATION                                   
!     COMMON SUBROUTINE CALL                                            
!c      CALL LTINV  (S,KDR,MJD)                                           
call ltinv  (s,kdr)
!     MATRIX S ARRANGEMENT 
if(kdr.eq.1) go to 260
!xx      DO 12 I=2,KDR
do 13 i=2,kdr
im1=i-1
do 12 j=1,im1
!xx   12 S(I,J)=S(J,I)
s(i,j)=s(j,i)
12 continue
13 continue
260 continue
!                                                                       
do  500 j=1,l
!     AV=A*V COMPUTATION:                                               
!                                                                       
!     +-                               -+     +-    -+                  
!     I     0     I       0   ...    0  I     I V(1) I                  
!     I     0     0       I   ...    0  I     I V(2) I                  
!     I     0     0       0   ...    0  I  *  I V(3) I                  
!     I     .     .       .   ...    .  I     I   .  I                  
!     I     .     .       .   ...    .  I     I   .  I                  
!     I -B(K) -B(K-1) -B(K-2) ... -B(1) I     I V(K) I                  
!     +-                               -+     +-    -+                  
!                                                                       
!       +-                                   -+                         
!       I  V(2)                               I                         
!     = I  V(3)                               I ,                       
!       I  V(4)                               I                         
!       I    .                                I                         
!       I    .                                I                         
!       I  -B(K)V(1)-B(K-1)V(2)-...-B(1)V(K)  I                         
!       +-                                   -+                         
!                                                                       
!     WHERE V(J)(J=1,K) ARE THE D-DIMENSIONAL SUB VECTORS OF V.         
!     A*V IS STORED IN THE VECTOR AV.                                   
!xx      DO  432 I1=1,D                                                    
!xx  432 W1(I1)=CST0
w1(1:d)=cst0
iz=k+1
kjz=0
do  436 jy=1,k
iz=iz-1
do  434 i1=1,d
do  433 i2=1,d
j2=kjz+i2
!xx  433 W1(I1)=W1(I1)+B(I1,I2,IZ)*V(J2)
w1(i1)=w1(i1)+b(i1,i2,iz)*v(j2)
433 continue
434 continue
kjz=kjz+d
436 continue
ii1=0
do  439 iy=1,k1
i1=ii1+1
i2=ii1+d
do  438 iz=i1,i2
j1=iz+d
av(iz)=v(j1)
438 continue
ii1=ii1+d
439 continue
do  441 i1=1,d
i2=kx+i1
av(i2)=w1(i1)
441 continue
!     CI=G*V; CONTROLLER INPUT                                          
!     COMMON SUBROUTINE CALL                                            
!c      CALL MULVER(G,V,CI,R,KD,MJD,MJKD)                                 
call mulver(g,v,ci,r,kd)
!     BCN=BC*CI                                                         
!     COMMON SUBROUTINE CALL                                            
!c      CALL MULVER(BC,CI,BCN,KD,R,MJKD,MJD)                              
call mulver(bc,ci,bcn,kd,r)
!     RANDOM VECTOR SI GENERATION                                       
do 302 i=1,kdr
bdmn=cstm6
do 303 jrn=1,12
!c  303 BDMN=BDMN+RANDOM(0)
!xx  303 BDMN=BDMN+RANDM(0,KK1,KK2,KK3,KK4)
bdmn=bdmn+randm(0,kk1,kk2,kk3,kk4)
303 continue
!xx  302 BDN(I)=BDMN
bdn(i)=bdmn
302 continue
!     RANDOM NORMAL VECTOR BDN GENERATION                               
!     COMMON SUBROUTINE CALL                                            
!c      CALL LTRVEC(S,BDN,SI,KDR,KDR,MJD,MJD)
call ltrvec(s,bdn,si,kdr,kdr)
!     COMMON SUBROUTINE CALL                                            
!c      CALL MULVER(BD,SI,BDN,KD,KDR,MJKD,MJD)
call mulver(bd,si,bdn,kd,kdr)
!     V=AV+BCN+BDN                                                      
do  309 i=1,kd
!xx  309 V(I) = AV(I)+BCN(I)+BDN(I)
v(i) = av(i)+bcn(i)+bdn(i)
309 continue
do  310 i=1,d
!xx  310 Y(I,J)=V(I)
y(i,j)=v(i)
310 continue
500 continue
!                                                                       
!     ************************************************                  
!     MEAN, VARIANCE AND STANDARD DEVIATION COMPUTATION                 
!     ************************************************                  
al=l
al=cst1/al
do  600 i=1,d
sum=cst0
do  501 j=1,l
!xx  501 SUM=SUM+Y(I,J)                                                    
sum=sum+y(i,j)
501 continue
!c      AV(I)=AL*SUM                                                      
avy(i)=al*sum
sum=cst0
do  502 j=1,l
!c  502 SUM=SUM+(Y(I,J)-AV(I))**2                                         
!xx  502 SUM=SUM+(Y(I,J)-AVY(I))**2
sum=sum+(y(i,j)-avy(i))**2
502 continue
si(i)=al*sum
s2(i)= dsqrt(si(i))
600 continue
!c      WRITE(6,902)                                                      
!c      WRITE(6,904)                                                      
!c      DO  610 I=1,D                                                     
!c  610 WRITE(6,905) I,AV(I),SI(I),S2(I)                                  
!                                                                       
!c      CALL FLCLS2(NFL)
!c  999 CONTINUE
return
!xx  800 FORMAT(6I5)
!xx  900 FORMAT('1  PROGRAM 74.3.2. OPTIMAL CONTROLLER DESIGN')            
!xx  901 FORMAT(1H ,'D (DIMENSION)=',I2,', K (ORDER)=',I3,', H (HORIZON)=',
!xx     A       I3,', L (LENGTH OF EXPERIMENTAL PERIOD)=',I5,              
!xx     A       ', R (DIMENSION OF CONTROL INPUT)=',I5)                    
!xx  902 FORMAT(//1H ,'OPTIMAL CONTROL SIMULATION')                        
!xx  903 FORMAT(////1H ,'COTROLLER GAIN G(',I3,',',I3,')')                 
!xx  904 FORMAT(1H ,' AVERAGE VALUE OF I-TH COMPONENT OF Y',9X,            
!xx     A'VARIANCE',20X,'STANDARD DEVIATION')                              
!xx  905 FORMAT(1H ,'      I = ',I4,D20.10,8X,D20.10,8X,D20.10)            
!xx  700 FORMAT(' MATRIX  B(',I3,')')                                   
!xx  701 FORMAT(' MATRIX  BX')                                          
!xx  702 FORMAT(' MATRIX  BC')                                          
!xx  703 FORMAT(' MATRIX  BD')                                          
!xx  704 FORMAT(' MATRIX  QX')                                          
!xx  706 FORMAT(' S: COVARIANCE MATRIX OF INNOVATION')                  
end
!                                                                       
!c      SUBROUTINE MULTRX(X,MX,NX,MJ1,MJ2,Y,MY,NY,MJ3,MJ4,                
!c     A                  Z,MZ,NZ,MJ5,MJ6,IS)                             
subroutine multrx(x,mx,nx,y,my,ny,z,mz,nz,is)
  use timsac_kinds, only: dp
  implicit none
!     MATRIX MULTIPLICATION                                             
!     X: MX*NX, ABSOLUTE DIMENSION MJ1*MJ2                              
!     Y: MY*NY, ABSOLUTE DIMENSION MJ3*MJ4                              
!     Z: MZ*NZ, ABSOLUTE DIMENSION MJ5*MJ6                              
!     IS=1: Z=X*Y                                                       
!     IS=2: Z=X'*Y                                                      
!     IS=3: Z=X*Y'                                                      
!xx      IMPLICIT REAL*8 (A-H,O-Z)                                         
!c      DIMENSION X(MJ1,MJ2),Y(MJ3,MJ4),Z(MJ5,MJ6)                        
!xx      DIMENSION X(MX,NX),Y(MY,NY),Z(MZ,NZ)
integer mx, nx, my, ny, mz, nz, is
real(dp) x(mx,nx), y(my,ny), z(mz,nz)
! local
integer i, j, k
real(dp) cst0
cst0=0.0d-00
if  (is.eq.2) go to 3050
if  (is.eq.3) go to 3100
!c      MZ=MX                                                             
!c      NZ=NY                                                             
do  3009 i=1,mx
do  3008 j=1,ny
z(i,j)=cst0
do  3007 k=1,nx
!xx 3007 Z(I,J)=Z(I,J)+X(I,K)*Y(K,J)
z(i,j)=z(i,j)+x(i,k)*y(k,j)
3007 continue
3008 continue
3009 continue
go to 3200
!c 3050 MZ=NX                                                             
!c      NZ=NY                                                             
3050 continue
do  3059 i=1,nx
do  3058 j=1,ny
z(i,j)=cst0
do  3057 k=1,mx
!xx 3057 Z(I,J)=Z(I,J)+X(K,I)*Y(K,J)
z(i,j)=z(i,j)+x(k,i)*y(k,j)
3057 continue
3058 continue
3059 continue
go to 3200
!c 3100 MZ=MX                                                             
!c      NZ=MY                                                             
3100 continue
do  3159 i=1,mx
do  3158 j=1,my
z(i,j)=cst0
do  3157 k=1,nx
!xx 3157 Z(I,J)=Z(I,J)+X(I,K)*Y(J,K)
z(i,j)=z(i,j)+x(i,k)*y(j,k)
3157 continue
3158 continue
3159 continue
3200 return
end
!                                                                       
!c      SUBROUTINE INVDET(X,XDET,MM,MJ)                                   
subroutine invdets(x,xdet,mm)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE                                                 
!     THIS SUBROUTINE COMPUTES THE INVERSE AND DETERMINANT OF           
!     UPPER LEFT MM X MM OF X.                                          
!     X: ORIGINAL MATRIX                                                
!     MM: DIMENSION OF UPPER LEFT OF X (SHOULD BE LESS THAN 11)         
!     XDET: DETERMINANT OF UPPER LEFT MM X MM OF X                      
!     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE                   
!     THE INVERSE MATRIX IS OVERWRITTEN ON THE ORIGINAL.                
!     NEXT STATEMENT SHOULD BE REPLACED BY                              
!     IMPLICIT COMPLEX*16(X)                                            
!     FOR COMPLEX VERSION.  ALSO STATEMENT NO.1 NEEDS MODIFICATION.     
!xx      IMPLICIT REAL*8(X)                                                
!c      DIMENSION X(MJ,MJ)                                                
!c      DIMENSION IDS(10)                                                 
!xx      DIMENSION X(MM,MM)
!xx      DIMENSION IDS(MM)
integer mm
real(dp) x(mm,mm), xdet
! local
integer i, j, jj, l, maxi, mm1, mmj, ids(mm)
real(dp) cst0, cst1, xmaxp, xc
cst0=0.0d-00
cst1=1.0d-00
xdet=cst1
do 10 l=1,mm
!     PIVOTING AT L-TH STAGE                                            
xmaxp=0.10000d-10
maxi=0
do 110 i=l,mm
!     FOR COMPLEX VERSION NEXT STATEMENT SHOULD BE REPLACED BY          
!     IF(CDABS(XMAXP).GE.CDABS(X(I,L))) GO TO 110                       
!xx    1 IF(DABS(XMAXP).GE.DABS(X(I,L))) GO TO 110                 
if(dabs(xmaxp).ge.dabs(x(i,l))) go to 110
xmaxp=x(i,l)
maxi=i
110 continue
ids(l)=maxi
if(maxi.eq.l) go to 120
if(maxi.gt.0) go to 121
xdet=cst0
go to 140
!     ROW INTERCHANGE                                                   
121 do 14 j=1,mm
xc=x(maxi,j)
x(maxi,j)=x(l,j)
!xx   14 X(L,J)=XC                                                         
x(l,j)=xc
14 continue
xdet=-xdet
120 xdet=cst1
xc=cst1/xmaxp
x(l,l)=cst1
do 11 j=1,mm
!xx   11 X(L,J)=X(L,J)*XC
x(l,j)=x(l,j)*xc
11 continue
do 12 i=1,mm
if(i.eq.l) go to 12
xc=x(i,l)
x(i,l)=cst0
do 13 j=1,mm
!xx   13 X(I,J)=X(I,J)-XC*X(L,J)
x(i,j)=x(i,j)-xc*x(l,j)
13 continue
12 continue
10 continue
if(mm.gt.1) go to 123
go to 140
!     COLUMN INTERCHANGE                                                
123 mm1=mm-1
do 130 j=1,mm1
mmj=mm-j
jj=ids(mmj)
if(jj.eq.mmj) go to 130
do 131 i=1,mm
xc=x(i,jj)
x(i,jj)=x(i,mmj)
!xx  131 X(I,MMJ)=XC
x(i,mmj)=xc
131 continue
130 continue
140 return
end
