! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine mulnosf(h,l,ip,sd,a,rs1,rs2,r)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.3.3   MULTIPLE UNOISE
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
!     THIS PROGRAM COMPUTES RELATIVE POWER CONTRIBUTIONS IN DIFFERENTIAL
!     AND INTEGRATED FORM, ASSUMING THE ORTHOGONALITY BETWEEN NOISE
!     SOURCES.
!     THE PROGRAM OPERATES ON THE OUTPUT OF PROGRAM 5.3.2 FPEC WITH
!     IL=0.
!     THE RESULTS ARE GIVEN AT FREQUIENCIES I/(2*H).
!
!xx      IMPLICIT REAL*8(A-H,O-W)
!xx      IMPLICIT COMPLEX*16(X-Z)
!xx      INTEGER H,H1
!      COMMON G,GR,GI,LG,H,JJF
!      DIMENSION SD(10,10),A(30,10,10),X(10,10)
!      DIMENSION G(31),RS(10,10),R(10,10)
!xx      DIMENSION SD(IP,IP),A(L,IP,IP),X(IP,IP)
!xx      DIMENSION G(L+1),RS1(IP,IP),RS2(IP,IP,H+1),R(IP,IP,H+1)
integer h, l, ip
real(dp) sd(ip,ip), a(l,ip,ip), rs1(ip,ip),&
&rs2(ip,ip,h+1), r(ip,ip,h+1)
! local
integer i, i1, ii, j, jf, jj, jjf, lg, h1
real(dp) g(l+1), gr, gi, cst0, cst1
complex(kind(0.0d0)) x(ip,ip), xdet
!     INPUT / OUTPUT DATA FILE OPEN
!	CALL SETWND
!	CALL FLOPN2(NFL)
!	IF (NFL.EQ.0) GO TO 999
!     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
!      MJ0=30
!      MJ1=10
cst0=0.0d-00
cst1=1.0d-00
!     H SPECIFICATION
!      READ(5,1) H
!     READING THE OUTPUT OF PROGRAM 5.3.2 FPEC WITH IL=0
!      READ(5,1) N,L,IP
!     SD INPUT
!      CALL REMATX(SD,IP,IP,1,MJ1,MJ1)
!     A INPUT
!      CALL REMAT3(A,L,IP,IP,1,MJ0,MJ1,MJ1)
!xx  310 H1=H+1
h1=h+1
!     SD NORMALIZATION
!xx      DO 100 I=1,IP
do 101 i=1,ip
do 100 j=1,ip
!  100 RS(I,J)=SD(I,J)/DSQRT(SD(I,I)*SD(J,J))
!xx  100 RS1(I,J)=SD(I,J)/DSQRT(SD(I,I)*SD(J,J))
rs1(i,j)=sd(i,j)/dsqrt(sd(i,i)*sd(j,j))
100 continue
101 continue
!     INITIAL CONDITION PRINT OUT
!      WRITE(6,59)
!      WRITE(6,60)
!      WRITE(6,61) H,N,L,IP
!      WRITE(6,161)
!      CALL SUBMPR(SD,IP,IP,MJ1,MJ1)
!     NORMALIZED SD PRINT OUT
!      WRITE(6,162)
!      CALL SUBMPR(RS,IP,IP,MJ1,MJ1)
!      CALL SUBMPR(RS,IP,IP,IP,IP)
!     A PRINT OUT
!      WRITE(6,420)
!      CALL PRMAT3(A,L,IP,IP,0,MJ0,MJ1,MJ1)
!xx  410 DO 10 JF=1,H1
do 10 jf=1,h1
jjf=jf
!     AF COMPUTATION
do 40 ii=1,ip
do 41 jj=1,ip
if(ii.ne.jj) go to 42
g(1)=cst1
go to 43
42 g(1)=cst0
43 do 45 i=1,l
i1=i+1
!xx   45 G(I1)=-A(I,II,JJ)
g(i1)=-a(i,ii,jj)
45 continue
!xx   44 LG=L
lg=l
!      CALL FGER1
call fger1(g,gr,gi,lg,h,jjf)
!xx      X(II,JJ)=DCMPLX(GR,GI)
x(ii,jj)=cmplx(gr,gi,kind(0.0d0))
41 continue
40 continue
!     INVERSE OF AF (COMPLEX) COMPUTATION
!      CALL INVDET(X,XDET,IP,MJ1)
!      CALL INVDET(X,XDET,IP)
call invdetc(x,xdet,ip)
!     RELATIVE POWER CONTRIBUTIONS COMPUTATION
!      CALL SUBNOS(X,SD,IP,RS,R,MJ1)
!      CALL SUBNOS(X,SD,IP,RS,R,IP)
call subnos(x,sd,ip,rs2(1,1,jf),r(1,1,jf),ip)
!     RELATIVE POWER CONTRIBUTIONS PRINT OUT
!      JFM1=JF-1
!      WRITE(6,65) JFM1
!      WRITE(6,165)
!      CALL SUBMPR(RS,IP,IP,10,10)
!      WRITE(6,166)
!      CALL SUBMPR(R,IP,IP,10,10)
10 continue
!	CALL FLCLS2(NFL)
!  999 CONTINUE
return
!    1 FORMAT(10I5)
!    2 FORMAT(4D20.10)
!   59 FORMAT(1H ,31HPROGRAM 5.3.3   MULTIPLE UNOISE)
!   60 FORMAT(1H ,17HINITIAL CONDITION)
!   61 FORMAT(1H ,2HH=,I5,5X,2HN=,I5,5X,2HL=,I5,5X,3HIP=,I5)
!  161 FORMAT(/1H ,7HSD(I,J))
!  162 FORMAT(/1H ,13HNORMALIZED SD)
!  420 FORMAT(/1H ,6HA(I,J))
!   65 FORMAT(///1H ,2HF=,I5)
!  165 FORMAT(/1H ,40HDIFFERENTIAL RELATIVE POWER CONTRIBUTION)
!  166 FORMAT(/1H ,38HINTEGRATED RELATIVE POWER CONTRIBUTION)
end subroutine
