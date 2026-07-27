! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine wnoisef(nra,ir,sd1,x2)
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM WNOISE
!     PROGRAM 5.5.3   WHITE NOISE GENERATOR
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
!     THIS PROGRAM GENERATES APPROXIMATELY GAUSSIAN VECTOR WHITE NOISE
!     TO BE USED AS INPUT W OF PROGRAM 5.5.2 OPTSIM.
!     ON TOP OF THE OUTPUT OF PROGRAM 5.3.2 FPEC, WHICH IS TO BE USED
!     AS INPUT TO PROGRAM 5.5.2 OPTSIM, ONE CARD WITH SPECIFICATION
!     OF THE LENGTH NRA OF WHITE NOISE RECORD TO BE GENERATED SHOULD
!     BE ADDED TO FORM THE INPUT TO THIS PROGRAM.
!     NRA: LENGTH OF WHITE NOISE RECORD TO BE GENERATED
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      PARAMETER (MJ0=100)
integer, parameter :: mj0=100
!c      REAL*4 RANDOM
!cc	 REAL*4	 RANDM
!xx      REAL*8 RANDM
!c      DIMENSION SD(5,5),A(10),Y(5),Z(5)
!c      DIMENSION X1(100,5)
!xx      DIMENSION SD1(IR,IR),SD(IR,IR),Y(IR),Z(IR)
!xx      DIMENSION X1(MJ0,IR),X2(IR,NRA)
integer nra, ir
real(dp) sd1(ir,ir), x2(ir,nra)
! local
integer i, i1, i2, i3, ii, iic, im1, ind, ind1, j, jj, k1, k2,&
&k3, k4
real(dp) randm, sd(ir,ir), y(ir), z(ir), x1(mj0,ir), cst0,&
&rc, rconst, xx, sum
!     INPUT / OUTPUT DATA FILE OPEN
!c      CHARACTER(100) DFNAM
!c      DFNAM = 'wnoise.out'
!c      CALL SETWND
!c      CALL FLOPN3(DFNAM,NFL)
!c      IF (NFL.EQ.0) GO TO 999
!     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
!c      MJ0=100
!c      MJ=5
cst0=0.0d-00
!     NRA SPECIFICATION
!c      READ(5,1) NRA
!     READING THE OUTPUTS OF PROGRAM 5.3.2 FPEC
!c      READ(5,1) N,M,IR,IL
!c      CALL REMATX(SD,IR,IR,1,MJ,MJ)
!xx      DO 5 I=1,IR
do 6 i=1,ir
do 5 j=1,ir
sd(i,j)=sd1(i,j)
5 continue
6 continue
!     FOLLOWING INPUT IS NONEFFECTIVE.
!c      IP=IR+IL
!c      MR=M*IR
!c      DO 8 JJ=1,MR
!c    8 READ(5,2) (A(II),II=1,IP)
!c      WRITE(6,60)
!c      WRITE(6,61)
!c      WRITE(6,62) NRA,N,M,IR
!c      WRITE(6,63)
!c      CALL SUBMPR(SD,IR,IR,MJ,MJ)
!c      WRITE(6,100)
!     MATRIX L COMPUTATION
!c      CALL LTINV(SD,IR,MJ)
call ltinv(sd,ir)
!     MATRIX L ARRANGEMENT
if(ir.eq.1) go to 260
!xx      DO 12 I=2,IR
do 13 i=2,ir
im1=i-1
do 12 j=1,im1
!xx   12 SD(I,J)=SD(J,I)
sd(i,j)=sd(j,i)
12 continue
13 continue
!     RANDOM NUMBER GENERATION
260 rc=4.0d-00/3.0d-00
rconst=dsqrt(rc)
!c      XX=RANDOM(1)
!cc      XX=RANDM(1)
xx=randm(1,k1,k2,k3,k4)
!c      IND=99
ind=mj0-1
ind1=ind+1
i1=0
i2=0
iic=0
160 i1=i2+1
i2=i1+ind
if(i2.le.nra) go to 130
i2=nra
130 do 14 i=i1,i2
ii=i-iic
do 15 j=1,ir
sum=cst0
do 20 jj=1,9
!c   20 SUM=SUM+RANDOM(0)
!cc   20 SUM=SUM+RANDM(0)
!xx   20 SUM=SUM+RANDM(0,K1,K2,K3,K4)
sum=sum+randm(0,k1,k2,k3,k4)
20 continue
sum=sum-4.5d-00
!xx   15 X1(II,J)=SUM*RCONST
x1(ii,j)=sum*rconst
15 continue
14 continue
!     WHITE NOISE GENERATION
do 16 i=i1,i2
ii=i-iic
do 17 j=1,ir
!xx   17 Y(J)=X1(II,J)
y(j)=x1(ii,j)
17 continue
!c      CALL LTRVEC(SD,Y,Z,IR,IR,MJ,MJ)
call ltrvec(sd,y,z,ir,ir)
do 18 j=1,ir
!xx   18 X1(II,J)=Z(J)
x1(ii,j)=z(j)
18 continue
16 continue
!     WHITE NOISE PRINT AND PUNCH OUT
!c      WRITE(6,101) IIC
i3=i2-iic
!c      CALL SUBMPR(X1,I3,IR,MJ0,MJ)
!xx      DO 40 I=I1,I2
do 41 i=i1,i2
ii=i-iic
!c   40 WRITE(7,3) (X1(II,J),J=1,IR)
do 40 j=1,ir
!xx   40 X2(J,I)=X1(II,J)
x2(j,i)=x1(ii,j)
40 continue
41 continue
iic=iic+ind1
if(i2.lt.nra) go to 160
!c      CALL FLCLS3(NFL)
!c  999 CONTINUE
return
!xx    1 FORMAT(10I5)
!xx    2 FORMAT(4D20.10)
!xx    3 FORMAT(6D12.3)
!xx   60 FORMAT(1H ,27HPROGRAM 5.5.3   WHITE NOISE)
!xx   61 FORMAT(1H ,17HINITIAL CONDITION)
!xx   62 FORMAT(1H ,4HNRA=,I5,5X,2HN=,I5,5X,2HM=,I5,5X,3HIR=,I5)
!xx   63 FORMAT(/1H ,7HSD(I,J))
!xx  100 FORMAT(////1H ,11HWHITE NOISE)
!xx  101 FORMAT(1H ,4HIIC=,I5,5X,11HX1(IIC+I,J))
end
