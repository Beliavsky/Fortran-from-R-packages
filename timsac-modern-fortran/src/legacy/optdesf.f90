! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine optdesf(ir,l,ns,m,q1,r,gr1,a,b,gi)
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM OPTDES
!     PROGRAM 5.5.1   OPTIMAL CONTROLLER DESIGN
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "DAINAMIKKU SISTEMU NO TOKEI-TEKI KAISEKI TO SEIGYO (STATISTICA
!         ANALYSIS AND CONTROL OF DYNAMIC SYSTEMS)" BY H. AKAIKE AND
!         T. NAKAGAWA, SAIENSU-SHA, TOKYO, 1972 (IN JAPANESE)
!-----------------------------------------------------------------------
!     THIS PROGRAM COMPUTES OPTIMAL CONTROLLER GAIN MATRIX FOR
!     A QUADRATIC CRITERION DEFINED BY TWO POSITIVE DEFINITE MATRICES
!     Q1 AND R.
!     THE OUTPUT OF FPEC COMPUTATION IS USED.
!     THE FIRST IR VARIABLES SHOULD BE CONTROLLED VARIABLES.
!     IR: NUMBER OF CONTROLLED VARIABLES
!     L: NUMBER OF MANIPULATED VARIABLES
!     NS: NUMBER OF D.P. STAGES
!     N: LENGTH OF ORIGINAL DATA
!     M: ORDER OF THE MODEL WHICH GIVES THE MINIMUM OF FPEC
!     MATRIX P: MI OR P
!     GI: GAIN
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION A(75,5),B(75,5),Q1(5,5),R(5,5)
!c      DIMENSION GI(5,75),GIT(75,5)
!c      DIMENSION GL(5,5),G3(5,75),GR(5,5),GLR(5,5)
!c      DIMENSION P(75,75),D(75)
!xx      DIMENSION A(M*IR,IR),B(M*IR,L),Q1(IR,IR),R(L,L)
!xx      DIMENSION GI(L,M*IR),GIT(M*IR,IR)
!xx      DIMENSION GL(L,L),G3(L,M*IR),GR(IR,IR),GLR(L,IR)
!xx      DIMENSION P(M*IR,M*IR),D(M*IR)
!xx      DIMENSION GR1(IR,IR)
integer ir, l, ns, m
real(dp) q1(ir,ir), r(l,l), gr1(ir,ir), a(m*ir,ir),&
&b(m*ir,l), gi(l,m*ir)
! local
integer i, i1, i2, iba, ii, ii0, iib, im1, ins, j, j1, j2, jj,&
&jjc, mm1, mm2, mr, mr1
real(dp) git(m*ir,ir), gl(l,l), g3(l,m*ir), gr(ir,ir),&
&glr(l,ir), p(m*ir,m*ir), d(m*ir), xdet, cst0
!     INPUT / OUTPUT DATA FILE OPEN
!c      CHARACTER(100) DFNAM
!c      DFNAM='optdes.out'
!c      CALL SETWND
!c      CALL FLOPN3(DFNAM,NFL)
!c      IF (NFL.EQ.0) GO TO 999
!     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
!c      MJ1=5
!c      MJ2=5
!c      MJ3=75
cst0=0.0d-00
!     INITIAL CONDITION INPUT AND OUTPUT
!c      READ(5,1) IR,L,NS
!c      CALL REMTSB(Q1,IR,1,MJ1)
!c      CALL REMTSB(R,L,1,MJ2)
!     READING THE OUTPUTS OF PROGRAM 5.3.2 FPEC
!c      READ(5,1) N,M
!c      CALL REMATX(GR,IR,IR,1,MJ1,MJ1)
!xx      DO 5 I=1,IR
do 6 i=1,ir
do 5 j=1,ir
gr(i,j)=gr1(i,j)
5 continue
6 continue
!     ABOVE INPUT IS NONEFFECTIVE.
!c      IM=0
!c      DO 8 JJ=1,M
!c      DO 9 I=1,IR
!c      II=IM+I
!c    9 READ(5,2) (A(II,J),J=1,IR),(B(II,J),J=1,L)
!c    8 IM=IM+IR
mr=m*ir
!c      WRITE(6,60)
!c      WRITE(6,61)
!c      WRITE(6,362) N,M
!c      WRITE(6,62) IR,L,NS
!c      WRITE(6,65)
!c      CALL SUBMPR(Q1,IR,IR,MJ1,MJ1)
!c      WRITE(6,66)
!c      CALL SUBMPR(R,L,L,MJ2,MJ2)
!c      WRITE(6,63)
!c      CALL SUBMPR(A,MR,IR,MJ3,MJ1)
!c      WRITE(6,64)
!c      CALL SUBMPR(B,MR,L,MJ3,MJ2)
!     CONSTANTS FOR COMPUTATION
mr1=mr-ir
mm1=m-1
mm2=m-2
iba=mm2*ir
!     INITIAL P COMPUTATION
!xx      DO 130 I=1,MR
do 132 i=1,mr
do 130 j=1,mr
!xx  130 P(I,J)=CST0
p(i,j)=cst0
130 continue
132 continue
!xx      DO 131 I=1,IR
do 133 i=1,ir
do 131 j=1,ir
!xx  131 P(I,J)=Q1(I,J)
p(i,j)=q1(i,j)
131 continue
133 continue
!
!     P ITERATION INS=1 TO NS
do 10 ins=1,ns
!     GI=(TRANSPOSE OF B)*P COMPUTATION
!c      CALL TRAMDL(B,P,GI,MR,L,MR,MJ3,MJ2,MJ3)
call tramdl(b,p,gi,mr,l,mr)
!     GL=GI*B COMPUTATION
!c      CALL MULTRB(GI,B,GL,L,MR,MJ2,MJ3)
call multrb(gi,b,gl,l,mr)
!xx      DO 20 I=1,L
do 21 i=1,l
do 20 j=1,i
!xx   20 GL(J,I)=GL(I,J)
gl(j,i)=gl(i,j)
20 continue
21 continue
!     GL=GL+R COMPUTATION
!c      CALL MATADL(GL,R,L,L,MJ2,MJ2)
call matadl(gl,r,l,l)
!     GL=(INVERSE OF GL) COMPUTATION
!c      CALL INVDET(GL,XDET,L,MJ2)
call invdet(gl,xdet,l,l)
!     G3=GL*GI COMPUTATION
!c      CALL MULPLY(GL,GI,G3,L,L,MR,MJ2,MJ2,MJ3)
call mulply(gl,gi,g3,l,l,mr)
!     D= DIAGONAL OF P
do 30 i=1,mr
!xx   30 D(I)=P(I,I)
d(i)=p(i,i)
30 continue
!     INTERMEDIATE MATRIX MI COMPUTATION
!     P=(TRANSPOSE OF GI)*G3 (LOWER TRIANGLE) COMPUTATION
!c      CALL MULTRL(GI,G3,P,L,MR,MJ2,MJ3)
call multrl(gi,g3,p,l,mr)
do 40 i=1,mr
!xx   40 P(I,I)=D(I)-P(I,I)
p(i,i)=d(i)-p(i,i)
40 continue
if(mr.eq.1) go to 260
!xx      DO 41 I=2,MR
do 42 i=2,mr
im1=i-1
do 41 j=1,im1
p(i,j)=p(j,i)-p(i,j)
!xx   41 P(J,I)=P(I,J)
p(j,i)=p(i,j)
41 continue
42 continue
!     GIT=MI*A COMPUTATION
!c  260 CALL MULPLY(P,A,GIT,MR,MR,IR,MJ3,MJ3,MJ1)
260 call mulply(p,a,git,mr,mr,ir)
!     GR=(TRANSPOSE OF A)*GIT COMPUTATION
!c      CALL MULTRL(A,GIT,GR,MR,IR,MJ3,MJ1)
call multrl(a,git,gr,mr,ir)
!xx      DO 150 I=1,IR
do 151 i=1,ir
do 150 j=1,i
!xx  150 GR(J,I)=GR(I,J)
gr(j,i)=gr(i,j)
150 continue
151 continue
!     GR=GR+Q1 COMPUTATION
!c      CALL MATADL(GR,Q1,IR,IR,MJ1,MJ1)
call matadl(gr,q1,ir,ir)
!     NEW P ARRANGEMENT
if(m.eq.1) go to 261
iib=iba
do 50 ii=1,mm1
ii0=m-ii
jjc=iib
do 51 jj=1,ii0
do 52 i=1,ir
i1=iib+i
i2=i1+ir
!X      DO 52 J=1,IR
do 53 j=1,ir
j1=jjc+j
j2=j1+ir
p(i2,j2)=p(i1,j1)
!xx   52 P(J2,I2)=P(I2,J2)
p(j2,i2)=p(i2,j2)
53 continue
52 continue
!xx   51 JJC=JJC-IR
jjc=jjc-ir
51 continue
!xx   50 IIB=IIB-IR
iib=iib-ir
50 continue
do 57 i=1,mr1
i1=i+ir
!xx      DO 57 J=1,IR
do 58 j=1,ir
p(i1,j)=git(i,j)
!xx   57 P(J,I1)=P(I1,J)
p(j,i1)=p(i1,j)
58 continue
57 continue
!xx  261 DO 56 I=1,IR
261 do 55 i=1,ir
do 56 j=1,i
p(i,j)=gr(i,j)
!xx   56 P(J,I)=P(I,J)
p(j,i)=p(i,j)
56 continue
55 continue
10 continue
!
!     GAIN COMPUTATION
!     GLR=-G3*A
!c      CALL MULPLY(G3,A,GLR,L,MR,IR,MJ2,MJ3,MJ1)
call mulply(g3,a,glr,l,mr,ir)
!xx      DO 110 I=1,L
do 112 i=1,l
do 110 j=1,ir
!xx  110 GI(I,J)=-GLR(I,J)
gi(i,j)=-glr(i,j)
110 continue
112 continue
if(m.eq.1) go to 262
!xx      DO 111 I=1,L
do 113 i=1,l
do 111 j=1,mr1
j1=ir+j
!xx  111 GI(I,J1)=-G3(I,J)
gi(i,j1)=-g3(i,j)
111 continue
113 continue
!     GAIN PRINT AND PUNCH OUT
262 continue
!c  262 WRITE(6,67)
!c      CALL SUBMPR(GI,L,MR,MJ2,MJ3)
!c      WRITE(7,1) N,M,IR,L,NS
!c      DO 232 I=1,IR
!c  232 WRITE(7,2) (Q1(I,J),J=1,IR)
!c      DO 233 I=1,L
!c  233 WRITE(7,2) (R(I,J),J=1,L)
!c      DO 230 I=1,MR
!c  230 WRITE(7,2) (A(I,J),J=1,IR)
!c      DO 231 I=1,MR
!c  231 WRITE(7,2) (B(I,J),J=1,L)
!c      DO 120 I=1,L
!c  120 WRITE(7,2) (GI(I,J),J=1,MR)
!c      CALL FLCLS3(NFL)
!c  999 CONTINUE
return
!xx    1 FORMAT(10I5)
!xx    2 FORMAT(4D20.10)
!xx   60 FORMAT(1H ,41HPROGRAM 5.5.1   OPTIMAL CONTROLLER DESIGN)
!xx   61 FORMAT(1H ,17HINITIAL CONDITION)
!xx   62 FORMAT(1H ,3HIR=,I5,5X,2HL=,I5,5X,3HNS=,I5)
!xx  362 FORMAT(1H ,2HN=,I5,5X,2HM=,I5)
!xx   63 FORMAT(//1H ,44HFIRST IR COLUMNS OF TRANSITION MATRIX (AI'S))
!xx   64 FORMAT(//1H ,19HGAMMA MATRIX (BI'S))
!xx   65 FORMAT(//1H ,10X,7HQ1(I,J))
!xx   66 FORMAT(//1H ,10X,6HR(I,J))
!xx   67 FORMAT(////1H ,13HGAIN MATRIX G)
end
!
!c	SUBROUTINE MULTRB(X,Y,Z,MM,NN,MJ1,MJ2)
subroutine multrb(x,y,z,mm,nn)
  use timsac_kinds, only: dp
  implicit none
!     Z=X*Y
!     Z: SYMMETRIC
!     (LOWER TRIANGLE OF UPPER LEFT MM X MM OF Z)=(UPPER LEFT MM X NN OF
!     X)*(UPPER LEFT NN X MM OF Y)
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!     (MJ2,MJ1): ABSOLUTE DIMENSION OF Y IN THE MAIN ROUTINE
!     (MJ1,MJ1): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION X(MJ1,MJ2),Y(MJ2,MJ1),Z(MJ1,MJ1)
!xx      DIMENSION X(MM,NN),Y(NN,MM),Z(NN,NN)
integer mm, nn
real(dp) x(mm,nn), y(nn,mm), z(nn,nn)
! local
integer i, j, k
real(dp) cst0, sum
cst0=0.0d-00
!xx      DO 10 I=1,MM
do 9 i=1,mm
do 10 j=1,i
sum=cst0
do 11 k=1,nn
!xx   11 SUM=SUM+X(I,K)*Y(K,J)
sum=sum+x(i,k)*y(k,j)
11 continue
!xx   10 Z(I,J)=SUM
z(i,j)=sum
10 continue
9 continue
return
end
!
!c      SUBROUTINE MULTRL(X,Y,Z,MM,NN,MJ1,MJ2)
subroutine multrl(x,y,z,mm,nn)
  use timsac_kinds, only: dp
  implicit none
!     TRANSPOSE MULTIPLY (LEFT)
!     Z=X'*Y
!     Z: SYMMETRIC
!     (LOWER TRIANGLE OF UPPER LEFT NN X NN OF Z)=(UPPER LEFT MM X NN OF
!     X)'*(UPPER LEFT MM X NN OF Y)
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X AND Y IN THE MAIN ROUTINE
!     (MJ2,MJ2): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION X(MJ1,MJ2),Y(MJ1,MJ2),Z(MJ2,MJ2)
!xx      DIMENSION X(MM,NN),Y(MM,NN),Z(NN,NN)
integer mm, nn
real(dp) x(mm,nn), y(mm,nn), z(nn,nn)
! local
integer i, j, k
real(dp) cst0, sum
cst0=0.0d-00
!xx      DO 10 I=1,NN
do 9 i=1,nn
do 10 j=1,i
sum=cst0
do 11 k=1,mm
!xx   11 SUM=SUM+X(K,I)*Y(K,J)
sum=sum+x(k,i)*y(k,j)
11 continue
!xx   10 Z(I,J)=SUM
z(i,j)=sum
10 continue
9 continue
return
end
