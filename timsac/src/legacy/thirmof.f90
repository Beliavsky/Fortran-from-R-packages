! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine thirmof(n,mh,xx,xmean,cc,cn,c)
  use timsac_kinds, only: dp
  implicit none
!
!c	PROGRAM THIRMO
!     PROGRAM 74.6.1.
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
!     THIS PROGRAM COMPUTES THE THIRD ORDER MOMENTS.
!     THE OUTPUT OF THIS PROGRAM M,MH,CC(I) AND C(I,J) ARE USED AS THE I
!     TO THE PROGRAM BISPEC (BI-SPECTRUM COMPUTATION).
!     THIS PROGRAM REQUIRES THE FOLLOWING INPUTS;
!     N: DATA LENGTH
!     MH: MAXIMUM LAG
!     DFORM: DATA FORMAT STATEMENT INFORMATION
!     X(I),I=1,N: ORIGINAL DATA
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      REAL*4 X,XMEAN
!c      DIMENSION X(2000),X1(2000),CC(51),C(51,51)
!c      DIMENSION CN(51),ITAB(7),CPR(7),IIND(7)
!c      DIMENSION DFORM(20)
!xx      DIMENSION X(N),X1(N),CC(MH+1),CN(MH+1),C(MH+1,MH+1)
!xx      DIMENSION XX(N)
integer n, mh
real(dp) xx(n), xmean, cc(mh+1), cn(mh+1), c(mh+1,mh+1)
! local
integer i, j, k, ki, kj, kl1, kl2, l1
real(dp) x(n), x1(n), cst0, cst1, an, bn, sum, t1, t2, cx0
!c	DATA IIND / 7*4HJ=   /
!
!     INPUT / OUTPUT DATA FILE OPEN
!c      CHARACTER(100) DFNAM
!c      CALL SETWND
!c      DFNAM='thirmo.out'
!c      CALL FLOPN3(DFNAM,NFL)
!c      IF (NFL.EQ.0) GO TO 999
!
!c      READ(5,1) N,MH
l1=mh+1
!     INPUT FORMAT SPECIFICATION
!c      READ(5,4) (DFORM(I),I=1,20)
!c    4 FORMAT(20A4)
!     ORIGINAL DATA INPUT AND PRINT OUT
!c      READ(5,DFORM) (X(I),I=1,N)
do 9 i=1,n
x(i)=xx(i)
9 continue
!c      WRITE(6,49)
!c      WRITE(6,50)
!c      WRITE(6,51)
!c      WRITE(6,52) N,MH
!c      WRITE(6,53)
!c      WRITE(6,55) ((X(I),I),I=1,N)
!
!     MEAN DELETION
cst0=0.0d-00
cst1 =1.0d-00
an=n
bn =cst1/an
sum =cst0
do 10 i=1,n
!xx   10 SUM=SUM+X(I)
sum=sum+x(i)
10 continue
xmean=bn*sum
do 11 i=1,n
!xx   11 X(I)=X(I)-XMEAN
x(i)=x(i)-xmean
11 continue
do 20 j=1,l1
kl1=n-j+1
t1=cst0
do 21 k=1,kl1
kj=k+j-1
!c      X1(K) =DBLE(X(KJ))*DBLE(X(K))
x1(k) =x(kj)*x(k)
!xx   21 T1=T1+X1(K)
t1=t1+x1(k)
21 continue
cc(j)=bn*t1
do 22 i=j,l1
kl2=n-i+1
t2 =cst0
do 23 k=1,kl2
ki=k+i-1
!xx   23 T2=T2+X(KI)*X1(K)
t2=t2+x(ki)*x1(k)
23 continue
c(i,j)=bn*t2
22 continue
20 continue
!     NORMALIZATION
cx0=cc(1)
!     COMMON SUBROUTINE CALL
call cornom(cc,cn,l1,cx0,cx0)
!     AUTO COVARIANCE PRINT OUT
!c      WRITE(6,1162)
!c      WRITE(6,162) N,MH,XMEAN
!c      WRITE(6,163)
!c      DO 360 I=1,L1
!c      IM1=I-1
!c  360 WRITE(6,155) IM1,CC(I),CN(I)
!     C(I,J) PRINT OUT
!c      WRITE(6,180)
!c      DO 300 I=1,L1
!c      JG=-1
!c      JS=0
!c      JC=0
!c  319 JC=JC+1
!c      IF  (JC.GT.I) GO TO 305
!c      JS=JS+1
!c      IF  (JS.LE.7) GO TO 310
!c      JS=7
!c  316 WRITE(6,190) (IIND(J),ITAB(J),J=1,JS)
!c  316 CONTINUE
!c      JG=JG+1
!c      IF  (JG) 317,318,317
!c  317 WRITE(6,192) (CPR(J),J=1,JS)
!c  317 CONTINUE
!c      GO TO 315
!c  318 IM1=I-1
!c      WRITE(6,191) IM1,(CPR(J),J=1,JS)
!c  315 IF  (JC.EQ.54321) GO TO 300
!c      JS=1
!c  310 ITAB(JS)=JC-1
!c      CPR(JS)=C(I,JC)
!c      GO TO 319
!c  305 JC=54321
!c      GO TO 316
!c  300 CONTINUE
!     AUTO COVARIANCE PUNCH OUT
!c      WRITE(7,1) N,MH
!c      WRITE(7,2) (CC(I),I=1,L1)
!     C(I,J) PUNCH OUT
!c      DO 320 I=1,L1
!c  320 WRITE(7,2) (C(I,J),J=1,I)
!c      CALL FLCLS3(NFL)
!c  999 CONTINUE
return
!xx    1 FORMAT(16I5)
!xx    2 FORMAT(4D20.10)
!xx   49 FORMAT(1H ,'PROGRAM 74.6.1. THIRMO')
!xx   50 FORMAT(1H ,'THIRD ORDER MOMENTS COMPUTATION')
!xx   51 FORMAT(/1H ,'INITIAL CONDITION')
!xx   52 FORMAT(/1H ,'N=',I5,5X,'MH=',I5)
!xx   53 FORMAT(/1H ,'ORIGINAL DATA')
!xx   55 FORMAT(5(F13.4,' (',I5,'), '))
!xx 1162 FORMAT(///1H ,'AUTOCOVARIANCE')
!c  162 FORMAT(1H ,'N=',I5,5X,'MH=',I5,5X,'MEAN=',E15.5)
!xx  162 FORMAT(/1H ,'N=',I5,5X,'MH=',I5,5X,'MEAN=',D15.5)
!xx  163 FORMAT(/1H ,4X,'I',6X,'CC(I)=E(X(T+I)X(T))',6X,'NORMALIZED')
!xx  180 FORMAT(///1H ,'C(I,J)=E(X(T+I)X(T+J)X(T))')
!xx  181 FORMAT(1H ,2I5,D17.5)
!xx  155 FORMAT(1H ,I5,D17.5,8X,D17.5)
!xx  190 FORMAT(/1H ,7X,7(6X,A2,I5,4X))
!xx  191 FORMAT(1H ,'I=',I5,7D17.5)
!xx  192 FORMAT(1H ,7X,7D17.5)
end
