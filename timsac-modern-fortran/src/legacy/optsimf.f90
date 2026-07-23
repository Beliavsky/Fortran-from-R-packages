! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine optsimf(ns,m,ir,l,a,b,g,w,x,y,xmean,ymean,xs2,ys2,&
&xs2mea,ys2mea,xvar,yvar)
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM OPTSIM
!     PROGRAM 5.5.2   OPTIMAL CONTROL SIMULATION
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
!     THIS PROGRAM PERFORMS OPTIMAL CONTROL SIMULATION FOR THE
!     CONTROLLER DESIGNED BY PROGRAM 5.5.1 AND EVALUATES THE MEANS AND
!     VARIANCES OF THE CONTROLLED AND MANIPULATED VARIABLES X AND Y.
!     FOLLOWING CONSTANTS SHOULD BE PROVIDED BESIDES THE OUTPUTS OF
!     PROGRAM 5.5.1   OPTIMAL CONTROLLER DESIGN TO START THIS PROGRAM.
!     NS: NUMBER OF STEPS OF SIMULATION
!     INTP=1: TO SUPPRESS HISTORY OUTPUT
!     INTP=2: TO PRINT OUT THE HISTORY
!     THE SEQUENCE OF NS IR-DIMENSIONAL VECTORS W, REPRESENTING WHITE
!     NOISE (OR IMPULSE) IS ALSO REQUIRED AS INPUT.
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	DIMENSION A(75,5),B(75,5),Q1(5,5),R(5,5),G(5,75)
!c	DIMENSION X(5),XS(5),XS2(5)
!c	DIMENSION Y(5),YS(5),YS2(5)
!c	DIMENSION W(5),Z(75),C(75)
!c	DIMENSION XMEAN(5),XS2MEA(5),XVAR(5)
!c	DIMENSION YMEAN(5),YS2MEA(5),YVAR(5)
!xx      DIMENSION A(IR*M,IR),B(IR*M,L),G(L,IR*M)
!xx      DIMENSION X(IR,NS),XS(IR),XS2(IR)
!xx      DIMENSION Y(L,NS),YS(L),YS2(L)
!xx      DIMENSION W(IR,NS),Z(IR*M),C(IR*M)
!xx      DIMENSION XMEAN(IR),XS2MEA(IR),XVAR(IR)
!xx      DIMENSION YMEAN(L),YS2MEA(L),YVAR(L)
integer ns, m, ir, l
real(dp) a(ir*m,ir), b(ir*m,l), g(l,ir*m), w(ir,ns),&
&x(ir,ns), y(l,ns), xmean(ir), ymean(l), xs2(ir),&
&ys2(l), xs2mea(ir), ys2mea(l), xvar(ir), yvar(l)
! local
integer i, ins, ipr, mr, mr1
real(dp) xs(ir), ys(l), z(ir*m), c(ir*m), cst0, cst1,&
&ans, bns
!     INPUT / OUTPUT DATA FILE OPEN
!c	CALL SETWND
!c	CALL FLOPN2(NFL)
!c	IF (NFL.EQ.0) GO TO 999
!     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
!c	MJ1=5
!c	MJ2=5
!c	MJ3=75
cst0=0.0d-00
!     INITIAL CONDITION INPUT AND OUTPUT
!c	READ(5,1) NS,INTP
!     READING THE OUTPUTS OF PROGRAM 5.5.1 OPTDES
!c	READ(5,1) N,M,IR,L
mr=m*ir
!c	CALL REMATX(Q1,IR,IR,1,MJ1,MJ1)
!c	CALL REMATX(R,L,L,1,MJ2,MJ2)
!c	CALL REMATX(A,MR,IR,1,MJ3,MJ1)
!c	CALL REMATX(B,MR,L,1,MJ3,MJ2)
!c	CALL REMATX(G,L,MR,1,MJ2,MJ3)
!c	WRITE(6,60)
!c	WRITE(6,61)
!c	WRITE(6,62) N,M,IR,L,NS
!c	WRITE(6,65)
!c	CALL SUBMPR(Q1,IR,IR,MJ1,MJ1)
!c	WRITE(6,66)
!c	CALL SUBMPR(R,L,L,MJ2,MJ2)
!c	WRITE(6,63)
!c	CALL SUBMPR(A,MR,IR,MJ3,MJ1)
!c	WRITE(6,64)
!c	CALL SUBMPR(B,MR,L,MJ3,MJ2)
!c	WRITE(6,67)
!c	CALL SUBMPR(G,L,MR,MJ2,MJ3)
!     INITIAL CONDITIONING
!xx      DO 6 I=1,IR
!c	X(I)=CST0
!xx      DO 5 J=1,NS
!xx    5 X(I,J)=CST0
!xx      XS(I)=CST0
!xx      XS2(I)=CST0
!xx    6 CONTINUE
!xx      DO 7 I=1,L
!c	Y(I)=CST0
!xx      DO 77 J=1,NS
!xx   77 Y(I,J)=CST0
!xx      YS(I)=CST0
!xx      YS2(I)=CST0
!xx    7 CONTINUE
!xx      DO 8 I=1,MR
!xx      C(I)=CST0
!xx    8 CONTINUE
x(1:ir,1:ns)=cst0
xs(1:ir)=cst0
xs2(1:ir)=cst0
y(1:l,1:ns)=cst0
ys(1:l)=cst0
ys2(1:l)=cst0
c(1:mr)=cst0
mr1=mr-ir
!     START OF SIMULATION
!     NOISE INPUT
do 10 ins=1,ns
!c	READ(5,3) (W(I),I=1,IR)
!     X COMPUTATION
!c	CALL VECADL(C,W,IR)
call vecadl(c,w(1,ins),ir)
do 9 i=1,ir
!c    9 X(I)=C(I)
!xx    9 X(I,INS)=C(I)
x(i,ins)=c(i)
9 continue
!     Y COMPUTATION
!c	CALL MULVER(G,C,Y,L,MR,MJ2,MJ3)
call mulver(g,c,y(1,ins),l,mr)
if(ins.eq.ns) go to 101
!c	CALL MULVER(A,X,Z,MR,IR,MJ3,MJ1)
call mulver(a,x(1,ins),z,mr,ir)
if(m.eq.1) go to 360
do 20 i=1,mr1
ipr=i+ir
!xx   20 Z(I)=Z(I)+C(IPR)
z(i)=z(i)+c(ipr)
20 continue
!c  360 CALL MULVER(B,Y,C,MR,L,MJ3,MJ2)
360 call mulver(b,y(1,ins),c,mr,l)
call vecadl(c,z,mr)
!     SUM AND SUM OF SQUARES COMPUTATION
!c  101 CALL VECADL(XS,X,IR)
!c	CALL VECADL(YS,Y,L)
101 call vecadl(xs,x(1,ins),ir)
call vecadl(ys,y(1,ins),l)
do 30 i=1,ir
!c   30 XS2(I)=XS2(I)+X(I)**2
!xx   30 XS2(I)=XS2(I)+X(I,INS)**2
xs2(i)=xs2(i)+x(i,ins)**2
30 continue
do 31 i=1,l
!c   31 YS2(I)=YS2(I)+Y(I)**2
!xx   31 YS2(I)=YS2(I)+Y(I,INS)**2
ys2(i)=ys2(i)+y(i,ins)**2
31 continue
!c	IF(INTP.EQ.1) GO TO 10
!     X,Y,W PRINT OUT
!c	WRITE(6,260) INS
!c	WRITE(6,261)
!c	CALL PRCOL2(X,W,1,IR,0)
!c	WRITE(6,262)
!c	CALL PRCOL1(Y,1,L,0)
10 continue
!     MEAN, MEAN SQUARE AND VARIANCE COMPUTATION
ans=ns
cst1=1.0d-00
bns=cst1/ans
do 40 i=1,ir
xmean(i)=bns*xs(i)
xs2mea(i)=bns*xs2(i)
!xx   40 XVAR(I)=XS2MEA(I)-XMEAN(I)**2
xvar(i)=xs2mea(i)-xmean(i)**2
40 continue
do 41 i=1,l
ymean(i)=bns*ys(i)
ys2mea(i)=bns*ys2(i)
!xx   41 YVAR(I)=YS2MEA(I)-YMEAN(I)**2
yvar(i)=ys2mea(i)-ymean(i)**2
41 continue
!c	WRITE(6,160)
!c	WRITE(6,161)
!c	CALL PRCOL4(XMEAN,XS2,XS2MEA,XVAR,1,IR,0)
!c	WRITE(6,163)
!c	WRITE(6,164)
!c	CALL PRCOL4(YMEAN,YS2,YS2MEA,YVAR,1,L,0)
!c	CALL FLCLS2(NFL)
!c  999 CONTINUE
return
!xx    1 FORMAT(10I5)
!xx    3 FORMAT(6D12.3)
!xx   60 FORMAT(1H ,42HPROGRAM 5.5.2   OPTIMAL CONTROL SIMULATION)
!xx   61 FORMAT(1H ,17HINITIAL CONDITION)
!xx   62 FORMAT(1H ,2HN=,I5,5X,2HM=,I5,5X,3HIR=,I5,5X,2HL=,I5,5X,3HNS=,I5)
!xx   63 FORMAT(//1H ,44HFIRST IR COLUMNS OF TRANSITION MATRIX (AI'S))
!xx   64 FORMAT(//1H ,19HGAMMA MATRIX (BI'S))
!xx   65 FORMAT(//1H ,7HQ1(I,J))
!xx   66 FORMAT(//1H ,6HR(I,J))
!xx   67 FORMAT(//1H ,13HGAIN MATRIX G)
!xx  160 FORMAT(//1H ,29X,4HX(I))
!xx  161 FORMAT(1H ,4X,1HI,5X,9HMEAN OF X,5X,11HSUM OF X**2,3X,12HMEAN OF X
!xx     A**2,2X,13HVARIANCE OF X)
!xx  163 FORMAT(//1H ,29X,4HY(I))
!xx  164 FORMAT(1H ,4X,1HI,5X,9HMEAN OF Y,5X,11HSUM OF Y**2,3X,12HMEAN OF Y
!xx     A**2,2X,13HVARIANCE OF Y)
!xx  260 FORMAT(/1H ,4HINS=,I5)
!xx  261 FORMAT(1H ,4X,1HI,5X,4HX(I),10X,4HW(I))
!xx  262 FORMAT(1H ,4X,1HI,5X,4HY(I))
end
!
subroutine vecadl(x,y,mm)
  use timsac_kinds, only: dp
  implicit none
!     X=X+Y (X,Y: VECTORS)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION X(MM),Y(MM)
integer mm
real(dp) x(mm), y(mm)
! local
integer i
do 10 i=1,mm
!xx   10 X(I)=X(I)+Y(I)
x(i)=x(i)+y(i)
10 continue
return
end
