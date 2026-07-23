! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine canarmf(n,lagh3,cyy,coef,ifpl1,sd,aic,oaic,mo,a,&
!x     *NC,MM1,MM2,V,Z,Y,XX,NDT,X3,X3MIN,MIN3,M1M,BETA,M1N,ALPHA,TMP,
!x     *MJ1,MJ2,IER)
&nc,mm1,mm2,v,z,y,xx,ndt,x3,x3min,min3,m1m,beta,m1n,alpha,mj1,mj2)
  use timsac_kinds, only: dp
  implicit none
!
!c	PROGRAM CANARM
!     PROGRAM 74.1.1. CANONICAL CORRELATION ANALYSIS OF SCALAR TIME SERI
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(1
!         BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
!         NO.5, MARCH 1975, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** FOR THE BASIC THEORY SEE "CANONICAL CORRELATION ANALYSIS OF TIM
!         AND THE USE OF AN INFORMATION CRITERION" BY H. AKAIKE, IN
!         "SYSTEM IDENTIFICATION: ADVANCES AND CASE STUDIES" R. K. MEHRA
!         D. G. LAINIOTIS EDS. ACADEMIC PRESS, NEW YORK, 1976
!-----------------------------------------------------------------------
!     THIS PROGRAM FITS AN AR-MA MODEL TO STATIONARY SCALAR TIME SERIES
!     THROUGH THE ANALYSIS OF CANONICAL CORRELATIONS
!     BETWEEN THE FUTURE AND PAST SETS OF OBSERVATIONS.
!     THE OUTPUTS OF THIS PROGRAM SHOULD BE ADDED TO THE INPUTS
!     TO THIS PROGRAM TO FORM AN INPUT TO THE PROGRAM AUTARM.
!
!     INPUTS REQUIRED:
!     (N,LAGH0): N, LENGTH OF ORIGINAL DATA Y(I) (I=1,N)
!                    LAGH0, MAXIMUM LAG OF COVARIANCE
!     CYY(I),I=0,LAGH0: AUTOCOVARIANCE SEQUENCE OF Y(I)
!
!     OUTPUTS:
!     NEWL: NEWL=1, FOR DIRECT INPUT TO PROGRAM AUTARM
!     M1M: ORDER OF AR
!     BETA(I)(I=1,M1M): AR-COEFFICIENTS
!     M1N: ORDER OF MA (=M1M-1)
!     ALPHA(I)(I=1,M1N): MA-COEFFICIENTS
!
!     THE AR-MA MODEL IS GIVEN BY
!     Y(N)+BETA(1)Y(N-1)+...+BETA(M1M)Y(N-M1M) = X(N)+ALPHA(1)X(N-1)+...
!                                                                      ...+ALPHA(M1N)X(N-M1N)
!
!c      PARAMETER (MJ1=50,MJ2=101)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION X(50),Y(50),Z(50)
!c      DIMENSION WL(50)
!c      DIMENSION XX(50),X3(50),BETA(50)
!c      DIMENSION CYY(1001)
!c      DIMENSION A(101),ALPHA(50)
!c      DIMENSION VC(101),VT(101)
!c      DIMENSION ST(101,50),T(101,50),V(50,50)
!c      COMMON /COM9/AST1
!c      DIMENSION AST1(5200)
!c      DIMENSION VV(50,50)
!
!xx      DIMENSION CYY(LAGH3)
!xx      DIMENSION COEF(MJ2)
!xx      DIMENSION SD(0:MJ1),AIC(0:MJ1),A(MJ1)
!xx      DIMENSION MM1(MJ1),MM2(MJ1),V(MJ1,MJ1,MJ1)
!xx      DIMENSION Z(MJ1,MJ1),Y(MJ1,MJ1),XX(MJ1,MJ1)
!xx      DIMENSION NDT(MJ1,MJ1),X3(MJ1,MJ1)
!xx      DIMENSION X3MIN(MJ1),MIN3(MJ1)
!xx      DIMENSION BETA(MJ1),ALPHA(MJ1)
!xx      DIMENSION WL(MJ1)
!xx      DIMENSION VC(MJ2),VT(MJ2)
!xx      DIMENSION ST(MJ2,MJ1),T(MJ2,MJ1)
!xx      DIMENSION AST1((MJ2-1)*MJ2/2)
!xx      DIMENSION VV(MJ1,MJ1)
integer n, lagh3, ifpl1, mo, nc, m1m, m1n, mj1, mj2, mm1(mj1),&
&mm2(mj1), ndt(mj1,mj1), min3(mj1)
real(dp) cyy(lagh3), coef(mj2), sd(0:mj1), aic(0:mj1),&
&oaic, a(mj1), v(mj1,mj1,mj1), z(mj1,mj1),&
&y(mj1,mj1),  xx(mj1,mj1), x3(mj1,mj1),&
&x3min(mj1), beta(mj1), alpha(mj1)
! local
integer i, ii, indx, j, j1, lagh0, m, m1, m2, m9, na, newl, ninew,&
&ninew0
real(dp) wl(mj1), vc(mj2), vt(mj2), st(mj2,mj1),&
&t(mj2,mj1), ast1((mj2-1)*mj2/2), vv(mj1,mj1),&
&cst0, cst1, cst2, cst9, an, em, en, andt, aii
!
!x      INTEGER*1 TMP(1)
!x      CHARACTER CNAME*80
!
!     INITIAL CLEARING
!c      DATA BETA /50*0.0D-00/, ALPHA/50*0.0D-00/
!c      DO 100 I=1,50
!c         BETA(I)=0.0D0
!c         ALPHA(I)=0.0D0
!c  100 CONTINUE
!
!     INPUT / OUTPUT DATA FILE OPEN
!c      CHARACTER(100) DFNAM
!c      CALL SETWND
!c      DFNAM='canarm.out'
!c      CALL FLOPN3(DFNAM,NFL)
!c      IF (NFL.EQ.0) GO TO 999
!c      IF (NFL.EQ.0) GO TO 999
!
!x      LU=3
!x      IER=0
!x      DO 100 I = 1,80
!x  100 CNAME(I:I) = ' '
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
!xcx            WRITE(*,*) ' ***  canarm temp FILE OPEN ERROR :',CNAME,IVAR
!x            IER=IVAR
!x            IFG=0
!x         END IF
!x      END IF
!
!     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
!c      MJ1=50
!c      MJ2=101
cst0=0.0d-00
cst1=1.0d-00
cst2=2.0d-00
cst9=9999.0d-00
!c      DO 18 I=1,1001
!c   18 CYY(I)=CST0
!     INITIAL CONDITION INPUT AND PRINT OUT
!c      READ(5,1) N,LAGH0
!c      LAGH3=LAGH0+1
lagh0=lagh3-1
!c      READ(5,2) (CYY(I),I=1,LAGH3)
!c      WRITE(6,11130)
!c      WRITE(6,11111) N,LAGH0
an=n
!c      IFPL=3.0D-00*DSQRT(AN)
!c      IFPL=MIN0(IFPL,49,LAGH0)
!c      IFPL=MIN0(IFPL,(MJ1-1),LAGH0)
!c      IFPL1=IFPL+1
!     INITIAL AUTOREGRESSIVE MODEL FITTING
!c      CALL NSICP(CYY,IFPL1,N,A,MO,SD,AIC)
na=(ifpl1*(ifpl1+1)/2)
call nsicp(cyy,lagh3,ifpl1,n,ast1,na,coef,&
!xx     *                 SD,AIC,A,MO,OAIC,IFG,LU)
&sd,aic,a,mo,oaic)
ninew=mo
m9=mo+1
!     NINEW0=0
indx=1
!     MATRIX MULTIPLICATION: ORTHO-NORMALIZATION OF VARIABLES BY USING
!     AR-MODELS OF SUCCESSIVELY INCREASING ORDER
!     VC=AST1*CYY COMPUTATION
do 210 i=1,m9
!xx  210 VC(I)=AST1(1)*CYY(I)
vc(i)=ast1(1)*cyy(i)
210 continue
!     VT=VC*AST1' COMPUTATION
!     AST1 OBTAINED BY SUBROUTINE NSICP
!c      CALL SVCMAT(VC,VT,M9)
call svcmat(vc,vt,m9,ast1,na)
do 220 i=1,m9
!xx  220 ST(I,1)=VT(I)
st(i,1)=vt(i)
220 continue
!     NINEW0: NUMBER OF ITERATIONS
ninew0=1
nc=0
500 m1=ninew0+1
nc=nc+1
m2=m9
m=m1+m2
!     VC=(M1-TH ROW OF AST1)*CYY
!c      CALL SVECT(CYY,VC,M9,M1,INDX)
!c      CALL SVCMAT(VC,VT,M9)
call svect(cyy,lagh3,ast1,na,vc,m9,m1,indx)
call svcmat(vc,vt,m9,ast1,na)
do 240 i=1,m9
!xx  240 ST(I,M1)=VT(I)
st(i,m1)=vt(i)
240 continue
!xx      DO 250 I=1,M9
do 251 i=1,m9
do 250 j=1,m1
!xx  250 T(I,J)=ST(I,J)
t(i,j)=st(i,j)
250 continue
251 continue
!     SVD OF T
!     T IS THE COVARIANCE MATRIX BETWEEN THE SETS OF THE
!     ORTHO-NORMALIZED FUTURE AND PAST VARIABLES
!     SINGULAR VALUES (Z) ARE THE CANONICAL CORRELATION COEFFICIENTS.
!     COMMON SUBROUTINE CALL
!c      CALL MSVD(T,VV,Z,M2,M1,MJ2,MJ1)
!c      CALL SVTR(VV,V,M1,MJ1)
call msvd(t,vv,z(1,nc),m2,m1,mj2,mj1)
call svtr(vv,v(1,1,nc),ast1,na,m1,mj1)
do 260 j=1,m1
!c  260 Y(J)=Z(J)*Z(J)
!xx  260 Y(J,NC)=Z(J,NC)*Z(J,NC)
y(j,nc)=z(j,nc)*z(j,nc)
260 continue
!     FUTURE CANONICAL WEIGHTS (V) PRINT OUT
!c      WRITE(6,6)
!c      WRITE(6,7) M1
!c      WRITE(6,8) M2
!c      WRITE(6,9) N
!c      WRITE(6,35)
!     COMMON SUBROUTINE CALL
!c      CALL SUBMPR(V,M1,M1,MJ1,MJ1)
mm1(nc)=m1
mm2(nc)=m2
!     TEST OF DEPENDENCE OF THE LAST PREDICTOR BY DIC (DIFFERENCE OF AIC
!     DIC(J) = AIC(J) - AIC(MAXIMUM J)
em=m
en=n
j=m1
wl(j+1)=cst1
!c   42 WL(J)=WL(J+1)*(CST1-Y(J))
42 wl(j)=wl(j+1)*(cst1-y(j,nc))
j=j-1
if(j.gt.0) go to 42
do 45 j=1,m1
if(wl(j).gt.cst0) go to 145
!c      X(J)=CST9
!c      XX(J)=CST9
xx(j,nc)=cst9
go to 45
!c  145 X(J)=-EN*DLOG(WL(J))
!c      XX(J)=-EN*DLOG(WL(J))
145 xx(j,nc)=-en*dlog(wl(j))
45 continue
!c      NDT=M1*M2
!c      ANDT=NDT
ndt(1,nc)=m1*m2
andt=ndt(1,nc)
!     DIC(J)=X3(J)
!c      X3(1)=XX(1)-CST2*ANDT
x3(1,nc)=xx(1,nc)-cst2*andt
!c      WRITE(6,49)
j=0
!c      WRITE(6,50) J,Z(1),Y(1),XX(1),NDT,X3(1)
!c      X3MIN=X3(1)
!c      MIN3=0
x3min(nc)=x3(1,nc)
min3(nc)=0
if(m1.lt.2) go to 4110
do 51 j=2,m1
j1=j-1
!c      NDT=(M1-J1)*(M2-J1)
!c      ANDT=NDT
!c      X3(J)=XX(J)-CST2*ANDT
ndt(j,nc)=(m1-j1)*(m2-j1)
andt=ndt(j,nc)
x3(j,nc)=xx(j,nc)-cst2*andt
!c   51 WRITE(6,50) J1,Z(J),Y(J),XX(J),NDT,X3(J)
51 continue
do 4300 j=2,m1
!     MINIMUM OF DIC SERCH
!c      IF(X3(J).GE.X3MIN) GO TO 4300
!c      X3MIN=X3(J)
!c      MIN3=J-1
if(x3(j,nc).ge.x3min(nc)) go to 4300
x3min(nc)=x3(j,nc)
min3(nc)=j-1
4300 continue
!c 4110 WRITE(6,4410) X3MIN,MIN3
4110 continue
!c      WRITE(6,11112)
if(ninew0.eq.ninew) go to 6999
!     DEPENDENCE ACCEPTED WHEN M1N-DIC IS NEGATIVE
!c      IF(X3(M1).GT.CST0) GO TO 110
if(x3(m1,nc).gt.cst0) go to 110
6999 m1m=m1-1
if(m1m.le.0) go to 110
!     BETA(AR-COEFF) COMPUTATION
!c      AII=CST1/V(M1,M1)
aii=cst1/v(m1,m1,nc)
do 5100 i=1,m1m
ii=m1-i
!c      BETA(II)=V(M1,I)*AII
beta(ii)=v(m1,i,nc)*aii
5100 continue
!     ALPHA(MA-COEFF) COMPUTATION
!c      CALL ALPHAS(A,M1M,BETA,ALPHA)
call alphas(coef,m1m,beta,alpha)
!     THE INPUTS TO THE PROGRAM AUTAMA PUNCH OUT
5400 newl=1
!c      WRITE(7,1) NEWL
!     BETA(AR-COEFF), ALPHA(MA-COEFF) PUNCH OUT
!c      WRITE(7,1) M1M
if  (m1m.le.0) go to 5200
!c      WRITE(7,2) (BETA(I),I=1,M1M)
5200 m1n=m1m-1
!c      WRITE(7,1) M1N
if  (m1n.le.0) go to 1100
!c      WRITE(7,2) (ALPHA(I),I=1,M1N)
go to 1100
110 if(ninew0.lt.ninew) go to 5300
m1m=0
go to 5400
5300 ninew0=ninew0+1
go to 500
!c 1100 CALL FLCLS3(NFL)
1100 continue
!x      IF (IFG.NE.0) CLOSE(LU)
return
!xx    1 FORMAT(16I5)
!xx    2 FORMAT(4D20.10)
!xx    6 FORMAT(//1H ,21HCANONICAL CORRELATION)
!xx    7 FORMAT(/1H ,'NUMBER OF PRESENT AND FUTURE VARIABLES',2X,
!xx     A'M1=',I5)
!xx    8 FORMAT(/1H ,'NUMBER OF PRESENT AND PAST VARIABLES',4X,
!xx     A'M2=',I5)
!xx    9 FORMAT(/1H ,'DATA LENGTH=N=',I6)
!xx   35 FORMAT(/1H ,'FUTURE SET CANONICAL WEIGHTS, ROWWISE')
!xx   49 FORMAT(/1H ,5X,8HORDER(P),4X,11HCANONICAL R,3X,9HR-SQUARED,3X,
!xx     A 10HCHI-SQUARE,3X,6HN.D.F.,2X,23HDIC (P)(=CHI**2-2*D.F.))
!xx   50 FORMAT(1H ,10X,I3,4X,F8.4,4X,F8.3,7X,F8.2,4X,I6,3X,F10.4)
!xx 4410 FORMAT(/1H ,'MINIMUM DIC(P) =',F12.2,1X,'ATTAINED AT P=',I5)
!xx 1130 FORMAT(/1H ,'STRUCTURAL CHARACTERISTIC VECTOR',
!xx     A' (H(I),I=1,P)')
!xx11111 FORMAT(/1H ,'INITIAL AUTO REGRESSIVE MODEL FITTING ',
!xx     A'BY THE MINIMUM AIC PROCEDURE. / N=',I5,',LAGH0=',I5)
!xx11112 FORMAT(1H ,5X,'THE VALUES OF CHI-SQUARE AND DIC (P) ',
!xx     A'CORRESPONDING TO CANONICAL R=','1.000 ',
!xx     A'SHOULD BE IGNORED')
!xx11130 FORMAT(/1H ,'PROGRAM 74.1.1. CANARM')
end
!
subroutine alphas(a,m1m,beta,alpha)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES ALPHA(MA-COEFFICIENTS).
!c      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION A(101)
!xx      DIMENSION A(M1M)
!xx      DIMENSION BETA(M1M),ALPHA(M1M)
integer m1m
real(dp) a(m1m), beta(m1m), alpha(m1m)
! local
integer i, ipm, k, km1, kmi
real(dp) sum
alpha(m1m)=0.0d-00
if  (m1m.le.1) go to 20
alpha(1)=beta(1)-a(1)
if(m1m.le.2) go to 20
ipm=m1m-1
do 10 k=2,ipm
km1=k-1
sum=0.0
do 11 i=1,km1
kmi=k-i
!xx   11 SUM=SUM-ALPHA(I)*A(KMI)
sum=sum-alpha(i)*a(kmi)
11 continue
alpha(k)=beta(k)-a(k)+sum
10 continue
!     ALPHA,BETA PRINT OUT
!c   20 WRITE(6,60) M1M
20 continue
!c      WRITE(6,61)
!c      DO  21 I=1,M1M
!c   21 WRITE(6,62) I,BETA(I),ALPHA(I)
return
!xx   60 FORMAT(//1H ,'M1M=',I5)
!xx   61 FORMAT(/1H ,4X,1HI,6X,17HBETA(I) (AR-COEF),4X,
!xx     A 18HALPHA(I) (MA-COEF))
!xx   62 FORMAT(1H ,I5,6X,F17.5,5X,F17.5)
end
!
!
!c      SUBROUTINE NSICP(CYY,L1,N,COEF,MO,OSD,OAIC)
subroutine nsicp(cyy,l3,l1,n,ast1,na,coef,sd,aic,aa,&
!xx     *MO,OAIC,IFG,LU)
&mo,oaic)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE
!     THIS SUBROUTINE FITS AUTOREGRESSIVE MODELS OF SUCCESSIVELY
!     INCREASING ORDER UP TO L(=L1-1).
!     INPUT:
!     CYY(I),I=0,L1; AUTOCOVARIANCE SEQUENCE
!     L1: L1=L+1, L IS THE UPPER LIMIT OF THE MODEL ORDER
!     N; LENGTH OF ORIGINAL DATA
!     OUT PUT:
!     COEF; AR-COEFFICIENTS
!     MO: ORDER OF AR
!     OSD: INNOVATION VARIANCE
!     OAIC: VALUE OF AIC
!     AST1: MATRIX OF AR-COEFFICIENTS (IN VECTOR FORM)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM9/AST1
!c      DIMENSION AST1(5200)
!c      DIMENSION A(101),B(101)
!xx      DIMENSION CYY(L3),COEF(L1)
!xx      DIMENSION AST1(NA)
!xx      DIMENSION A(L1),B(L1)
!xx      DIMENSION AIC(0:L1),SD(0:L1),AA(L1)
integer l3, l1, n, na, mo
real(dp) cyy(l3), ast1(na), coef(l1), sd(0:l1), aic(0:l1),&
&aa(l1), oaic
! local
integer i, ian, ian1, ian2, ii, im, inx, jj, jj0, jjl, jjl1, l,&
&lan1, lan2, lm, m, mp1, nfc
real(dp) a(l1), b(l1), cst0, cst1, cst2, cst20, cst01,&
&cst05, am, an, anfc, ran, scalh, se, sdr, d, d2,&
&const, dlsd
!c      REAL*4 AX,BL,STA,DASH,PLUS
!c      REAL*4 FFFF
!c      REAL*4  F(41) / 41*1H  /, AMES(41) / 41*1H- /
!xx      CHARACTER FFFF
!xx      CHARACTER F(41),AMES(41)
!c      DATA AX,BL,STA,DASH,PLUS/1H!,1H ,1H*,1H-,1H+/
!xx      DATA F,AMES/ 41*' ', 41*'-' /
cst0=0.0d-00
cst1=1.0d-00
cst2=2.0d-00
cst20=20.0d-00
cst05=0.05d-00
cst01=0.00001d-00
l=l1-1
!c      SD=CYY(1)
sd(1)=cyy(1)
inx=1
!c      AST1(1)=CST1/DSQRT(SD)
ast1(1)=cst1/dsqrt(sd(1))
an=n
!c      OAIC=AN*DLOG(SD)
!c      OSD=SD
oaic=an*dlog(sd(1))
mo=0
!c      SD0=OSD
!c      AIC0=OAIC
sd(0)=sd(1)
aic(0)=oaic
!     INITIAL CONDITION PRINT OUT
!c  991 WRITE(6,1100)
!xx  991 CONTINUE
ran=cst1/dsqrt(an)
scalh=cst20
!xx      JJ0=SCALH+CST1
!xx      JJL=SCALH*CST2+CST1
jj0=int(scalh+cst1)
jjl=int(scalh*cst2+cst1)
jjl1=jjl-1
!xx      AMES(1)=PLUS
!xx      AMES(11)=PLUS
!xx      AMES(JJ0)=PLUS
!xx      AMES(JJ0+10)=PLUS
!xx      AMES(JJL)=PLUS
!xx      IAN=SCALH*(RAN+CST05)
ian=int(scalh*(ran+cst05))
ian1=ian+jj0
ian2=2*ian+jj0
lan1=-ian+jj0
lan2=-2*ian+jj0
!c	WRITE(6,26100)
!c	WRITE(6,26101)
!c	WRITE(6,261)
!c	WRITE(6,26102)
!c	WRITE(6,262)
!c	WRITE(6,264) (AMES(J),J=1,JJL)
!x      IF (IFG.NE.0) THEN
!x	 WRITE(LU,261)
!x	 WRITE(LU,262)
!x	 WRITE(LU,264) (AMES(J),J=1,JJL)
!x      END IF
!c	WRITE(6,859) MO,OSD,OAIC
!xx      F(JJ0)=AX
!xx      F(IAN1)=AX
!xx      F(IAN2)=AX
!xx      F(LAN1)=AX
!xx      F(LAN2)=AX
!c	WRITE(6,861) (F(J),J=1,JJL)
!x      IF (IFG.NE.0) WRITE(LU,264) (F(J),J=1,JJL)
se=cyy(2)
!     ITERATION START
do 400 m=1,l
!c	SDR=SD/CYY(1)
sdr=sd(m)/cyy(1)
if(sdr.ge.cst01) go to 399
!c	WRITE(6,2600)
go to 402
399 mp1=m+1
!c	D=SE/SD
d=se/sd(m)
a(m)=d
d2=d*d
!c	SD=(CST1-D2)*SD
!c	CONST=CST1/DSQRT(SD)
sd(m)=(cst1-d2)*sd(m)
const=cst1/dsqrt(sd(m))
am=m
!c	AIC=AN*DLOG(SD)+CST2*AM
aic(m)=an*dlog(sd(m))+cst2*am
!
!
!c	DLSD=DLOG(SD)
dlsd=dlog(sd(m))
if(m.eq.1) go to 410
!     A(I) COMPUTATION
lm=m-1
do 420 i=1,lm
a(i)=a(i)-d*b(i)
420 continue
410 do 460 i=1,m
ii=m+1-i
inx=inx+1
!xx  460 AST1(INX)=-A(II)*CONST
ast1(inx)=-a(ii)*const
460 continue
inx=inx+1
ast1(inx)=const
do 421 i=1,m
im=mp1-i
!xx  421 B(I)=A(IM)
b(i)=a(im)
421 continue
!     M,SD,AIC	PRINT OUT
if(a(m).lt.cst0) go to 300
!xx      NFC=SCALH*(A(M)+CST05)
nfc=int(scalh*(a(m)+cst05))
go to 310
!xx  300 NFC=SCALH*(A(M)-CST05)
300 nfc=int(scalh*(a(m)-cst05))
310 anfc=nfc
!xx      JJ=ANFC+SCALH+CST1
jj=int(anfc+scalh+cst1)
!xx      FFFF=F(JJ)
!xx      F(JJ)=STA
!c	WRITE(6,860) M,SD,AIC,A(M)
aa(m)=a(m)
!c	WRITE(6,861) (F(J),J=1,JJL)
!x      IF (IFG.NE.0) WRITE(LU,264) (F(J),J=1,JJL)
!xx      F(JJ)=FFFF
!c  990 IF(OAIC.LT.AIC) GO TO 440
!c	OAIC=AIC
!xx  990 IF(OAIC.LT.AIC(M)) GO TO 440
if(oaic.lt.aic(m)) go to 440
oaic=aic(m)
!c	OSD=SD
mo=m
do 430 i=1,m
!xx  430 COEF(I)=-A(I)
coef(i)=-a(i)
430 continue
440 if(m.eq.l) go to 400
se=cyy(m+2)
do 441 i=1,m
!xx  441 SE=SE-B(I)*CYY(I+1)
se=se-b(i)*cyy(i+1)
441 continue
sd(m+1)=sd(m)
400 continue
402 continue
!     MO, COEF(I) OUT PUT
!c	WRITE(6,870) OAIC,MO
!c	WRITE(6,1871)
!c	WRITE(6,871)
!c	CALL SUBVCP(COEF,MO)
!c  699 F(JJ0)=BL
!xx      F(IAN1)=BL
!xx      F(IAN2)=BL
!xx      F(LAN1)=BL
!xx      F(LAN2)=BL
!xx      AMES(JJ0)=DASH
!xx      AMES(JJ0+10)=DASH
!xx      AMES(JJL)=DASH
return
!xx26100 FORMAT(/1H ,16X,'SD(M)',15X,'AIC(M)',13X,'A(M)')
!xx26101 FORMAT(1H ,4X,'M',11X,'INNOVATION',10X,'AIC(M)=    ',8X,
!xx     A'PARTIAL AUTO-')
!c  261 FORMAT(1H ,16X,'   VARIANCE',9X,'N*DLOG(SD(M))+2*M',
!c     A      '  CORRELATION    ',
!c     A      'PARTIAL CORRELATION (LINES SHOW +SD AND +2SD)')
!xx  261 FORMAT(' PARTIAL CORRELATION (LINES SHOW +/-SD AND +/-2SD)')
!xx26102 FORMAT(1H ,102X,'_',7X,'_')
!c  262 FORMAT(1H ,69X,'-1',19X,'0',19X,'1')
!c  264 FORMAT(1H ,70X,41A1)
!xx  262 FORMAT(' -1',19X,'0',19X,'1')
!xx  264 FORMAT(2X,41A1)
!xx  859 FORMAT(1H ,I5,2X,2D20.5)
!xx  860 FORMAT(1H ,I5,2X,3D20.5)
!xx  861 FORMAT(1H ,70X,41A1)
!xx  960 FORMAT(/1H ,5X,'A(I)')
!xx  870 FORMAT(/1H ,'MINIMUM AIC(M)=',D12.5,2X,'ATTAINED AT M=',I5)
!xx  980 FORMAT(/1H ,5X,'COEF(I)')
!xx 1100 FORMAT(/1H ,'AIC(M)=N*DLOG(SD)+2*M')
!xx  871 FORMAT(/1H ,'AR-COEFFICIENTS')
!xx 1871 FORMAT(/1H ,'AR MODEL: Y(N)+AR(1)Y(N-1)+...+AR(M)Y(N-M)=X(N)')
!xx 2600 FORMAT(/1H ,'ACCURACY OF COMPUTATION LOST')
end
!
!c	SUBROUTINE SVCMAT(VC,VT,M9)
subroutine svcmat(vc,vt,m9,ast1,na)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES VT=VC*AST1'.
!     AST1 IS AN OUTPUT OF NSICP.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM9/AST1
!c      DIMENSION AST1(5200)
!xx      DIMENSION AST1(NA)
!xx      DIMENSION VC(M9),VT(M9)
integer m9, na
real(dp) vc(m9), vt(m9), ast1(na)
! local
integer i, inx, k
real(dp) cst0, sum
cst0=0.0d-00
inx=0
do  10 i=1,m9
sum=cst0
do 11 k=1,i
inx=inx+1
!xx   11 SUM=SUM+VC(K)*AST1(INX)
sum=sum+vc(k)*ast1(inx)
11 continue
vt(i)=sum
10 continue
return
end
!
subroutine svtr(vv,v,ast1,na,m1,mj1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES FUTURE SET CANONICAL WEIGHTS DEFINED BY
!     V=VV'*AST1.
!     AST1 IS AN OUTPUT OF NSICP.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM9/AST1
!c      DIMENSION AST1(5200),ISUM1(50)
!xx      DIMENSION VV(MJ1,MJ1)
!xx      DIMENSION AST1(NA),ISUM1(M1)
!xx      DIMENSION V(MJ1,MJ1)
integer na, m1, mj1
real(dp) vv(mj1,mj1), v(mj1,mj1), ast1(na)
! local
integer i, ijk, inx, isum, isum1(m1), j, k, ll
real(dp) cst0, sum
cst0=0.0d-00
isum=0
do 15 i=1,m1
isum=isum+i
!xx   15 ISUM1(I)=ISUM
isum1(i)=isum
15 continue
do 10 i=1,m1
do 11 j=1,m1
sum=cst0
ll=isum1(j)
inx=0
do 12 k=j,m1
ijk=ll+inx
sum=sum+vv(k,i)*ast1(ijk)
inx=inx+k
12 continue
v(i,j)=sum
11 continue
10 continue
return
end
!
subroutine svect(cyy,l3,ast1,na,vc,m9,m1,indx)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES VC=(M1-TH ROW OF AST1)*(CYY MATRIX)
!     AST1 IS AN OUTPUT OF NSICP.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM9/AST1
!c      DIMENSION AST1(5200)
!c      DIMENSION CYY(1001),VC(M9)
!xx      DIMENSION AST1(NA)
!xx      DIMENSION CYY(L3),VC(M9)
integer l3, na, m9, m1, indx
real(dp) cyy(l3), ast1(na), vc(m9), cst0
! local
integer i, ii, is, ism
cst0=0.0d-00
!xx      DO  10 I=1,M9
!xx   10 VC(I)=CST0
vc(1:m9)=cst0
do 20 is=1,m1
indx=indx+1
ism=is-1
do  30 i=1,m9
ii=ism+i
vc(i)=vc(i)+ast1(indx)*cyy(ii)
30 continue
20 continue
return
end
