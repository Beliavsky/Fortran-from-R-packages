! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine nonstf(n,istp,data0,nm,lagh,jp0,coef0,va0,aic0,daic21,&
&daic,k01,kount2,sxx)
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM NONST
!     PROGRAM 74.4.1. NON-STATIONARY POWER SPECTRUM ANALYSIS
!-----------------------------------------------------------------------
!     ** DESIGNED BY T. OZAKI, THE INSTITUTE OF STATISTICAL MATHEMATICS,
!     ** PROGRAMMED BY T. OZAKI, THE INSTITUTE OF STATISTICAL MATHEMATIC
!         TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(2
!         BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
!         NO.6 MARCH 1976, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** FOR THE BASIC THEORY SEE "ON THE FITTING OF NON-STATIONARY
!         AUTOREGRESSIVE MODELS IN TIME SERIES ANALYSIS" BY T. OZAKI AND
!         IN "PROC. 8TH HAWAII INTERNATIONAL CONFERENCE ON SYSTEM SCIENCE
!         WESTERN PERIODICALS NORTH HOLLYWOOD, CALIF., 1975
!-----------------------------------------------------------------------
!     THIS PROGRAM LOCALLY FITS AUTOREGRESSIVE MODELS TO NON-STATIONARY
!     TIME SERIES BY AIC CRITERION.
!     POWER SPECTRA FOR STATIONARY SPANS ARE GRAPHICALLY PRINTED OUT.
!     THE FOLLOWING INPUTS ARE REQUIRED;
!         N: LENGTH OF DATA
!         ISTP : LENGTH OF THE BASIC LOCAL SPAN
!         DFORM : INPUT FORMAT SPECIFICATION IN ONE CARD, FOR EXAMPLE,'(8
!         (X(I),I=1,N) : ORIGINAL DATA.
!
!     NSG : NUMBER OF SEGMENTS OF FREQUENCY AXIS +1
!xx      PARAMETER (NSG=121)
integer, parameter :: nsg=121
!
!xx      IMPLICIT REAL*8 (A-H,O-Z)
!c      DIMENSION ACV0(101),ACV1(101),ACV2(101),COEF0(50),COEF1(50)
!c      DIMENSION COEF2(50)
!xx      DIMENSION ACV0(LAGH+1),ACV1(LAGH+1),ACV2(LAGH+1)
!xx      DIMENSION COEF0(LAGH,NM),COEF1(LAGH),COEF2(LAGH)
!c      DIMENSION DATA1(9000),SXX(121)
!c      REAL*4 DATA(9000)
!c      REAL*4 DATA1(9000)
!c      REAL*4 DFORM(20)
!c      REAL*4 SXX(121)
!xx      DIMENSION DATA0(N),DATA1(N),SXX(NSG,NM)
!xx      INTEGER IDX(5)
!xx      DIMENSION CN(LAGH+1)
!xx      DIMENSION JP0(NM)
!xx      DIMENSION VA0(NM),AIC0(NM),DAIC21(NM),DAIC(NM)
!xx      DIMENSION K01(NM),KOUNT2(NM)
integer n, istp, nm, lagh, jp0(nm), k01(nm), kount2(nm)
real(dp) data0(n), coef0(lagh,nm), va0(nm), aic0(nm),&
&daic21(nm), daic(nm), sxx(nsg,nm)
! local
integer i, ii, ik0, ik1, ik5, ip0, ip1, ip2, j, jj, kount0,&
&kount1, lagh1, ml1, ni, nj, nml, nsg1, idx(5)
real(dp) acv0(lagh+1), acv1(lagh+1), acv2(lagh+1),&
&coef1(lagh), coef2(lagh), data1(n), cn(lagh+1),&
&aistp, aip0, aip1, aip2, cst2, cst4, ani, anj,&
&stp, va1, va2, aic1, aic2, xmean, z
!
!x      INTEGER*1  TMP(1)
!x      CHARACTER  CNAME*80
!
!     INPUT / OUTPUT DATA FILE OPEN
!c      CALL SETWND
!c      CALL FLOPN2(NFL)
!c      IF (NFL.EQ.0) GO TO 9999
!
!x      IER=0
!x      LU=3
!x      DO 10 I = 1,80
!x   10 CNAME(I:I) = ' '
!x      I = 1
!x      IFG = 1
!x      DO WHILE( (IFG.EQ.1) .AND. (I.LE.80) )
!x         IF ( TMP(I).NE.ICHAR(' ') ) THEN
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
!xcx            WRITE(*,*) ' ***  nonst temp FILE OPEN ERROR :',CNAME,IVAR
!x            IER=IVAR
!x            IFG=0
!x         END IF
!x      END IF
!
!     INITIALIZATION
!c      LAGH=50
lagh1=lagh+1
kount0=0
kount1=0
!c      KOUNT2=0
daic21(1)=0
daic(1)=0
!c      ML=49
!     ML : HIGHEST ALLOWABLE ORDER OF AUTOREGRESSIVE MODEL
ml1 = lagh
!     INITIAL  CONDITION  INPUT
!c      READ(5,1) N,ISTP
!     ORIGINAL	DATA  INPUT
!c      READ(5,4) (DFORM(I),I=1,20)
!c      READ(5,DFORM) (DATA(I),I=1,N)
!c      KOUNT2=ISTP
!c      KOUNT1=KOUNT2
kount2(1)=istp
kount1=kount2(1)
do 50 i=1,istp
!xx50    DATA1(I)=DATA0(I)
data1(i)=data0(i)
50 continue
!     COMON  SUBROUTINE	 CALL
!c	CALL AUTCOR(DATA1,ACV0,ISTP,LAGH)
!     COMON  SUBROUTINE	 CALL
!c	CALL ICPOUT(ACV0,COEF0,VA0,ML,ISTP,IP0,LAGH)
call autcorf(data1,istp,acv0,cn,lagh1,xmean)
call sicp(acv0,coef0,va0(1),ml1,istp,ip0)
aistp=istp
aip0=ip0
cst2=2.0d-00
!c	AIC0=AISTP*(DLOG(VA0))+CST2*(AIP0+CST2)
aic0(1)=aistp*(dlog(va0(1)))+cst2*(aip0+cst2)
ni=istp
!c      WRITE(6,8010)
!c      WRITE(6,8001) NI
!c      WRITE(6,8002) IP0
!c      WRITE(6,8003)
jp0(1)=ip0
!x      IF (IFG.NE.0) WRITE(LU,8001) NI
!x      IF (IFG.NE.0) WRITE(LU,8002) JP0(1)
idx(1)=1
do 60 i=1,4
j=i+1
!xx60    IDX(J)=IDX(I)+1
idx(j)=idx(i)+1
60 continue
70 ik1=idx(1)
ik5=idx(5)
ik0=ik5-ip0
if(ik0.ge.0) go to 90
!c      WRITE(6,8004) (IDX(I),I=1,5),(COEF0(I),I=IK1,IK5)
do 80 i=1,5
!xx80    IDX(I)=IDX(I)+5
idx(i)=idx(i)+5
80 continue
go to 70
90 ik0=5-ik0
!c      WRITE(6,8006) (IDX(I),I=1,IK0)
!c      WRITE(6,8005) (COEF0(I),I=IK1,IP0)
!c95    WRITE(6,8007) VA0
!c      WRITE(6,8008) AIC0
!c      K01=KOUNT0+1
k01(1)=kount0+1
!c      WRITE(6,8009) K01,KOUNT2
!x      IF (IFG.NE.0) WRITE(LU,8009) K01(1),KOUNT2(1)
!     COMON  SUBROUTINE	 CALL
!c      CALL NRASPE(VA0,COEF0,Z,SXX,IP0,0,ISTP,120,KOUNT0,KOUNT2)
!c      CALL SPEGRH(SXX)
nsg1=nsg-1
call nraspe(va0(1),coef0(:,1),[z],ip0,0,nsg1,sxx)
!x      IF (IFG.NE.0) CALL SPEGRH(LU,SXX,NSG)
nml=1
100 j=ni+1
jj=j+istp
!c      KOUNT2=KOUNT2+ISTP
!c      IF(KOUNT2.GT.N) GOTO 999
if((kount2(nml)+istp).gt.n) goto 999
kount2(nml+1)=kount2(nml)+istp
do 200 i=1,istp
ii=kount1+i
!xx200   DATA1(I)=DATA0(II)
data1(i)=data0(ii)
200 continue
nj=ni+istp
anj=nj
ani=ni
stp=istp
!     COMON  SUBROUTINE	 CALL
!c      CALL AUTCOR(DATA1,ACV2,ISTP,LAGH)
!     COMON  SUBROUTINE	 CALL
!c      CALL ICPOUT(ACV2,COEF2,VA2,ML,ISTP,IP2,LAGH)
call autcorf(data1,istp,acv2,cn,lagh1,xmean)
call sicp(acv2,coef2,va2,ml1,istp,ip2)
aip0=ip0
aip2=ip2
cst4=4.0d-00
!c      AIC2=ANI*DLOG(VA0)+STP*DLOG(VA2)+CST2*(AIP0+AIP2+CST4)
aic2=ani*dlog(va0(nml))+stp*dlog(va2)+cst2*(aip0+aip2+cst4)
do 500 j=1,nj
jj=kount0+j
!xx500   DATA1(J)=DATA0(JJ)
data1(j)=data0(jj)
500 continue
!c      CALL AUTCOR(DATA1,ACV1,NJ,LAGH)
!c      CALL ICPOUT(ACV1,COEF1,VA1,ML,NJ,IP1,LAGH)
call autcorf(data1,nj,acv1,cn,lagh1,xmean)
call sicp(acv1,coef1,va1,ml1,nj,ip1)
aip1=ip1
aic1=anj*dlog(va1)+cst2*(aip1+cst2)
!c      DAIC=AIC2-AIC1
nml=nml+1
daic(nml)=aic2-aic1
if(aic2.lt.aic1) goto 800
ni=nj
ip0=ip1
!c      VA0=VA1
!c      AIC0=AIC1
va0(nml)=va1
aic0(nml)=aic1
do 700 i=1,ip0
!c700	COEF0(I)=COEF1(I)
!xx700   COEF0(I,NML)=COEF1(I)
coef0(i,nml)=coef1(i)
700 continue
do 750 i=1,lagh
!xx750   ACV0(I)=ACV1(I)
acv0(i)=acv1(i)
750 continue
!c      WRITE(6,8010)
!c      WRITE(6,8001) NI
!c      WRITE(6,8002) IP0
!c      WRITE(6,8003)
jp0(nml)=ip0
!x      IF (IFG.NE.0) WRITE(LU,8001) NI
!x      IF (IFG.NE.0) WRITE(LU,8001) JP0(NML)
idx(1)=1
do 760 i=1,4
j=i+1
!xx760   IDX(J)=IDX(I)+1
idx(j)=idx(i)+1
760 continue
770 ik1=idx(1)
ik5=idx(5)
ik0=ik5-ip0
if(ik0.ge.0) go to 790
!c      WRITE(6,8004) (IDX(I),I=1,5),(COEF0(I),I=IK1,IK5)
do 780 i=1,5
!xx780   IDX(I)=IDX(I)+5
idx(i)=idx(i)+5
780 continue
go to 770
790 ik0=5-ik0
!c      WRITE(6,8006) (IDX(I),I=1,IK0)
!c      WRITE(6,8005) (COEF0(I),I=IK1,IP0)
!c795	WRITE(6,8007) VA0
!c      WRITE(6,8008) AIC0
!c      WRITE(6,8011) DAIC
daic21(nml)=daic(nml)
daic(nml)=daic(nml)/anj
!c      WRITE(6,8013) DAIC
!c      KOUNT1=KOUNT2
!c      K01=KOUNT0+1
kount1=kount2(nml)
k01(nml)=kount0+1
!c      WRITE(6,8009) K01,KOUNT2
!     COMON  SUBROUTINE	 CALL
!c      CALL NRASPE(VA0,COEF0,Z,SXX,IP0,0,NI,120,KOUNT0,KOUNT2)
!c      CALL SPEGRH(SXX)
call nraspe(va0(nml),coef0(:,nml),[z],ip0,0,nsg1,sxx(:,nml))
!x      IF (IFG.NE.0) CALL SPEGRH(LU,SXX(1,NML),NSG)
goto 100
800 ni=istp
ip0=ip2
!c      VA0=VA2
!c      AIC0=AIC2
va0(nml)=va2
aic0(nml)=aic2
do 1000 i=1,ip0
!c1000        COEF0(I)=COEF2(I)
!xx1000  COEF0(I,NML)=COEF2(I)
coef0(i,nml)=coef2(i)
1000 continue
do 1050 i=1,lagh
!xx1050  ACV0(I)=ACV2(I)
acv0(i)=acv2(i)
1050 continue
!c      WRITE(6,8010)
!c      WRITE(6,8000)
!c      WRITE(6,8001) NI
!c      WRITE(6,8002) IP0
!c      WRITE(6,8003)
jp0(nml)=ip0
!x      IF (IFG.NE.0) WRITE(LU,8001) NI
!x      IF (IFG.NE.0) WRITE(LU,8002) JP0(NML)
idx(1)=1
do 1060 i=1,4
j=i+1
!xx1060  IDX(J)=IDX(I)+1
idx(j)=idx(i)+1
1060 continue
1070 ik1=idx(1)
ik5=idx(5)
ik0=ik5-ip0
if(ik0.ge.0) go to 1090
!c      WRITE(6,8004) (IDX(I),I=1,5),(COEF0(I),I=IK1,IK5)
do 1080 i=1,5
!xx1080  IDX(I)=IDX(I)+5
idx(i)=idx(i)+5
1080 continue
go to 1070
1090 ik0=5-ik0
!c      WRITE(6,8006) (IDX(I),I=1,IK0)
!c      WRITE(6,8005) (COEF0(I),I=IK1,IP0)
!c1095      WRITE(6,8007) VA0
!c      WRITE(6,8008) AIC0
!c      WRITE(6,8011) DAIC
daic21(nml)=daic(nml)
daic(nml)=daic(nml)/anj
!c      WRITE(6,8013) DAIC
kount0=kount1
!c      KOUNT1=KOUNT2
!c      K01=KOUNT0+1
kount1=kount2(nml)
k01(nml)=kount0+1
!c      WRITE(6,8009) K01,KOUNT2
!     COMON  SUBROUTINE	 CALL
!c      CALL NRASPE(VA0,COEF0,Z,SXX,IP0,0,NI,120,KOUNT0,KOUNT2)
!c      CALL SPEGRH(SXX)
call nraspe(va0(nml),coef0(:,nml),[z],ip0,0,nsg1,sxx(:,nml))
!x      IF (IFG.NE.0) CALL SPEGRH(LU,SXX(1,NML),NSG)
goto 100
999 continue
!x      IF (IFG.NE.0) CLOSE(LU)
!c      CALL FLCLS2(NFL)
!c9999      CONTINUE
!xx1     FORMAT(20I5)
!xx4     FORMAT(20A4)
!xx8000  FORMAT(1H ,'SWITCHED TO A NEW MODEL HERE ]')
!xx8001  FORMAT(//1H ,'DATA  LENGTH  FOR  CURRENT  MODEL =',I6)
!xx8002  FORMAT(//,1H ,7X,'ORDER  OF  CURRENT  MODEL =',I6)
!xx8003  FORMAT(//,1H ,'AUTOREGRESSIVE  COEFFICIENTS')
!xx8004  FORMAT(1H ,15X,I3,22X,I3,22X,I3,22X,I3,22X,I3,/,1H ,5D25.10)
!xx8005  FORMAT(1H ,5D25.10)
!xx8006  FORMAT(1H ,15X,I3,22X,I3,22X,I3,22X,I3,22X,I3)
!xx8007  FORMAT(//,1H ,12X,'INNOVATION  VARIANCE =',D18.10)
!xx8008  FORMAT(//,1H ,29X,'AIC =',D18.10)
!xx8009  FORMAT(//,1H ,'THIS  MODEL  IS  FITTED  TO  THE  DATA  FROM',I6,
!xx     C'-TH  POINT  TO',I6,'-TH  POINT')
!xx8010  FORMAT(/////////////////////)
!xx8011  FORMAT(//1H ,21X,'AIC2 - AIC1 =',D18.10)
!xx8013  FORMAT(//,1H ,26X,'DAIC/N=',D18.10)
return
end
!
!
subroutine sicp(cxx,coef,osd,l1,n,mo)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE FITS AUTOREGRESSIVE MODELS
!     X(N)=A(1)X(N-1)+...+A(M)X(N-M)+E(N)
!     OF SUCCESSIVELY INCREASING ORDER UP TO L(=L1-1).
!     INPUT:
!     CXX(I),I=0,L1; AUTOCOVARIANCE SEQUENCE
!     L1; L1=L+1, L IS THE UPPER LIMIT OF THE MODEL ORDER
!     N; LENGTH OF ORIGINAL DATA
!     OUT PUT:
!     MO; ORDER OF AR
!     OSD; INNOVATION VARIANCE
!     COEF; AR-COEFFICIENTS
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION CXX(L1),COEF(L1)
!c	DIMENSION A(500),B(500)
integer l1, n, mo
real(dp) cxx(l1), coef(l1), osd
! local
integer i, im, l, lm, m, mp1
real(dp) a(l1-1),b(l1-1), cst1, cst2, sd, an, oaic, se,&
&d, d2, am, dlsd, aic
cst1=1.0d-00
cst2=2.0d-00
l=l1-1
sd=cxx(1)
an=n
oaic=an*dlog(sd)
osd=sd
mo=0
se=cxx(2)
!     ITERATION START
do 400 m=1,l
mp1=m+1
d=se/sd
a(m)=d
d2=d*d
sd=(cst1-d2)*sd
am=m
aic=an*dlog(sd)+cst2*am
dlsd=dlog(sd)
if(m.eq.1) go to 410
!
!     A(I) COMPUTATION
lm=m-1
do 420 i=1,lm
!xx  420 A(I)=A(I)-D*B(I)
a(i)=a(i)-d*b(i)
420 continue
410 do 421 i=1,m
im=mp1-i
!xx  421 B(I)=A(IM)
b(i)=a(im)
421 continue
!
!     M,SD,FICP PRINT OUT
!     WRITE(6,860) M
!     WRITE(6,861) SD,FICP
!     WRITE(6,960)
!     DO 961 I=1,M
! 961 WRITE(6,662) I,A(I)
!
if(oaic.lt.aic) go to 440
oaic=aic
osd=sd
mo=m
do 430 i=1,m
!xx  430 COEF(I)=A(I)
coef(i)=a(i)
430 continue
440 if(m.eq.l) go to 400
se=cxx(m+2)
do 441 i=1,m
!xx  441 SE=SE-B(I)*CXX(I+1)
se=se-b(i)*cxx(i+1)
441 continue
400 continue
!xx  699 RETURN
return
!xx    1 FORMAT(16I5)
!xx    2 FORMAT(4D20.10)
!xx  658 FORMAT(1H ,'INITIAL CONDITION')
!xx  659 FORMAT(1H ,2HN=,I5,5X,2HL=,I5)
!xx  660 FORMAT(1H ,'AUTO COVARIANCE')
!xx  661 FORMAT(1H ,4X,1HI,11X,6HCXX(I))
!xx  662 FORMAT(1H ,I5,D17.5)
!xx  757 FORMAT(1H ,2HM=,4X,1H )
!xx  758 FORMAT(1H ,'I.C.(P)=',D17.5)
!xx  860 FORMAT(1H ,2HM=,I5)
!xx  861 FORMAT(1H ,3HSD=,D17.5,5X,8HI.C.(P)=,D17.5)
!xx  960 FORMAT(1H ,4X,1HI,13X,4HA(I))
!xx870   FORMAT(1H ,'MINIMUM I.C.(P)=',D12.5,2X,'ATTAINED AT M=',I5)
!xx 1100 FORMAT(/1H ,'I.C.(P)=N*DLOG(SD)+2*M')
end
