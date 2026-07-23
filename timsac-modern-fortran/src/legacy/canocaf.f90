! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine canocaf(ir,inw,n,lagh1,ip0,ccv,l,aic,oaic,mo,osd,aao,&
&nc,n1,n2,vv,z,y,xx,ndt,x3,x3min,min3,f,m1nh,nh,g,iaw,vf,&
&lmax,mj0,mj1)
  use timsac_kinds, only: dp
  implicit none
!
!c	PROGRAM CANOCA
!     PROGRAM 74.2.1. CANONICAL CORRELATION ANALYSIS OF VECTOR TIME SERI
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
!     THIS PROGRAM DOES CANONICAL CORRELATION ANALYSIS OF AN IR-DIMENSIO
!     MULTIVARIATE TIME SERIES Y(I) (I=1,N).
!
!     FIRST AR-MODEL IS FITTED BY THE MINIMUM  A I C  PROCEDURE.
!     THE RESULTS ARE USED TO ORTHO-NORMALIZE THE PRESENT AND PAST VARIA
!     THE PRESENT AND FUTURE VARIABLES ARE TESTED SUCCESSIVELY TO DECIDE
!     ON THE DEPENDENCE OF THEIR PREDICTORS. WHEN THE LAST DIC (AN INFOR
!     CRITERION) IS NEGATIVE THE PREDICTOR OF THE VARIABLE IS DECIDED
!     TO BE LINEARLY DEPENDENT ON THE ANTECEDENTS.
!
!     THE STRUCTURAL CHARACTERISTIC VECTOR H OF THE CANONICAL MARKOVIAN
!     REPRESENTATION AND THE ESTIMATE OF THE TRANSITION MATRIX F, IN
!     VECTOR FORM, ARE PUNCHED OUT. THE ESTIMATE OF THE INPUT MATRIX G A
!     THE COVARIANCE MATRIX C OF THE INNOVATION, OBTAINED BY USING
!     THE F-MATRIX AND THE AR-MODEL, ARE ALSO PUNCHED OUT.
!
!     INPUTS REQUIRED:
!     IR:		 DIMENSION OF Y(I)
!     INW(K)(K=1,IP):	 INW(K)=J MEANS THAT THE K-TH COMPONENT OF Y(I)
!          IS THE J-TH COMPONENT OF THE ORIGINAL RECORD Z(I) USED FOR
!          THE COMPUTATION OF THE COVARIANCE SEQUENCE CZZ(I).
!
!     ALSO THE FOLLOWING INPUTS ARE REQUESTED BY THE SUBROUTINE	 R E C O
!     THESE INPUTS CAN BE OBTAINED AS THE OUTPUTS OF THE TIMSAC
!     PROGRAM 5.1.2 ( MULCOR ).
!     N:		 DATA LENGTH
!     LAGH:		 MAXIMUM LAG OF COVARIANCE
!     IP0:		 DIMENSION OF THE ORIGINAL RECORD
!     CZZ(I)(I=0,LAGH):	 COVARIANCE MATRIX SEQUENCE OF Z(I) GIVEN IN THE
!         OF THE SUCCESSION OF THE SEQUENCES OF COVARIANCES BETWEEN THE
!         C-TH COMPONENTS OF Z(I), ARRANGED IN THE ORDER (R=1,C=1),(R=2,
!         (R=2,C=1),(R=1,C=2),(R=3,C=3),(R=3,C=1),(R=1,C=3),(R=3,C=2),(R
!         ...... EACH SEQUENCE HAS THE HEADING (R,C).
!         THE (J+1)ST ELEMENT OF THE COVARIANCE SEQUENCE WITH R=K AND C=
!         IS AN ESTIMATE OF E(Z(I+J,K)*Z(I,L)), WHERE Z(I,K) DENOTES
!         THE K-TH COMPONENT OF Z(I).
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION Y(91),Z(91),WL(91)
!c      DIMENSION CV(25,7,7)
!c      DIMENSION XX(91),X3(91),F(91),VF(637)
!c      DIMENSION U(91,46),V(46,46),VV(46,46)
!c      DIMENSION AST1(91,7,7)
!c      DIMENSION NH(91),IH(91)
!c      DIMENSION FL(46,46)
!c      DIMENSION RGT(91,91),FRG(46,91)
!c      DIMENSION VTG(91),SR(46),ZZ(46),SF(46),SFRG(91)
!c      DIMENSION INW(7),C1(7,7)
!c      DIMENSION OSD(7,7),AO(13,7,7)
!c      DIMENSION G(46,7)
!xx      DIMENSION Y(MJ1,MJ1),Z(MJ1,MJ1),WL(MJ1)
!xx      DIMENSION CV(LMAX*2+1,IR,IR)
!xx      DIMENSION XX(MJ1,MJ1),X3(MJ1,MJ1),F(MJ1,MJ1),VF(MJ1*MJ1)
!xx      DIMENSION U(MJ1,MJ1),V(MJ1,MJ1),VV(MJ1,MJ1,MJ1)
!xx      DIMENSION AST1((LMAX+1)*LMAX,IR,IR)
!xx      DIMENSION NH(MJ1),IH(MJ1)
!xx      DIMENSION FL(MJ1,MJ1)
!xx      DIMENSION RGT(MJ1,MJ1),FRG(MJ1,MJ1)
!xx      DIMENSION VTG(MJ1),SR(MJ1),ZZ(MJ1),SF(MJ1),SFRG(MJ1)
!xx      DIMENSION INW(IR),C1(IP0,IP0)
!xx      DIMENSION OSD(IR,IR),AO(MJ0,IR,IR),AAO(MJ0,IR,IR)
!xx      DIMENSION G(MJ1,IR)
!
!xx      DIMENSION CCV(LAGH1,IP0,IP0)
!xx      DIMENSION AIC(MJ0)
!xx      DIMENSION N1(MJ1),N2(MJ1),NDT(MJ1,MJ1)
!xx      DIMENSION X3MIN(MJ0*IR),MIN3(MJ0*IR)
integer ir, n, lagh1, ip0, l, mo, nc, m1nh, iaw, lmax, mj0, mj1,&
&inw(ir), n1(mj1), n2(mj1), ndt(mj1,mj1), min3(mj0*ir),&
&nh(mj1)
real(dp) ccv(lagh1,ip0,ip0), aic(mj0), oaic, osd(ir,ir),&
&aao(mj0,ir,ir), vv(mj1,mj1,mj1), z(mj1,mj1),&
&y(mj1,mj1), xx(mj1,mj1), x3(mj1,mj1),&
&x3min(mj0*ir), f(mj1,mj1), g(mj1,ir), vf(mj1*mj1)
! local
integer i, id, idb, ii, ijdb, ijdb1, il, ip, ires, isw, j, j1,&
&jdb, jj, jres, k, l1, lcv, lcv1, lmax2, m, m1, m1m, m2,&
&mm1, mmmh, mp1, ninew, ninew0, ninewi, ih(mj1)
real(dp) wl(mj1), cv(lmax*2+1,ir,ir), u(mj1,mj1),&
&v(mj1,mj1), ast1((lmax+1)*lmax,ir,ir),&
&fl(mj1,mj1), rgt(mj1,mj1), frg(mj1,mj1),&
&vtg(mj1), sr(mj1), zz(mj1), sf(mj1), sfrg(mj1),&
&c1(ip0,ip0), ao(mj0,ir,ir), cst0, cst1, cst2, xl,&
&yl, r11c, fct, ss, fsrinp, ssfr2, sg, em, en,&
&andt, aii
!
!c      COMMON /COM9/AST1
!c      COMMON /COM10/CV
!c      COMMON /COM11/RGT
!c      COMMON /COM19/U
!c      COMMON /COM20/V
!c      COMMON /COM21/VV
!c      COMMON /COM25/FL
!c      EQUIVALENCE (AST1(1,1,1),VF(1))
!
!     INPUT / OUTPUT DATA FILE OPEN
!c      CHARACTER(100) DFNAM
!c      CALL SETWND
!c      DFNAM='canoca.out'
!c      CALL FLOPN3(DFNAM,NFL)
!c      IF (NFL.EQ.0) GO TO 999
!
!c      MJ1=46
!c      MJ2=91
!c      MJ=7
!c      MJ0=13
nc = 0
cst0=0.0d-00
cst1=1.0d-00
cst2=2.0d-00
!c      DO 3111 I=1,46
!c      DO 3112 J=1,7
do 3111 i=1,mj1
do 3112 j=1,ir
g(i,j)=cst0
3112 continue
!---
do 3113 j=1,mj1
f(i,j)=cst0
3113 continue
!---
3111 continue
!     INITIAL CONDITION INPUT
!c      READ(5,1) IR
il=0
ip=ir
!c      READ(5,1) (INW(I),I=1,IP)
!c      WRITE(6,11114)
!c      WRITE(6,11111)
!c      WRITE(6,11113) IR
!c      LCV=24
lcv=lmax*2
lcv1=lcv+1
!     AUTO COVARIANCE INPUT
!c      CALL RECOVA(LCV1,IP0,N)
!xx      DO 90 I=1,MIN(LCV1,LAGH1)
!xx         DO 90 J=1,IP
do 92 i=1,min(lcv1,lagh1)
do 91 j=1,ip
do 90 k=1,ip
cv(i,j,k) = ccv(i,j,k)
90 continue
91 continue
92 continue
!     L (UPPER BOUND OF AR-ORDER) COMPUTATION
!c      LMAX=12
xl=2*n
xl=dsqrt(xl)
yl=ir
!xx      J=XL/YL
j=int(xl/yl)
l= min0(lmax,j)
l1=l+1
!     MATRIX ARRANGEMENT
do 10 ii=1,lcv1
!xx      DO 20 I=1,IP0
do 11 i=1,ip0
do 20 j=1,ip0
!xx   20 C1(I,J)=CV(II,I,J)
c1(i,j)=cv(ii,i,j)
20 continue
11 continue
!     MATRIX REARRANGEMENT BY INW
!     COMMON SUBROUTINE CALL
!c      CALL REARRA(C1,INW,IP0,IP,MJ)
call rearra(c1,inw,ip0,ip)
!xx      DO 21 I=1,IP
do 12 i=1,ip
do 21 j=1,ip
!xx   21 CV(II,I,J)=C1(I,J)
cv(ii,i,j)=c1(i,j)
21 continue
12 continue
10 continue
!c      WRITE(6,259) (INW(I),I=1,IP)
iaw=0
id=ir
lmax2=(lmax+1)*lmax
!     AR-MODEL FITTING BY THE MINIMUM AIC PROCEDURE
!     COMMON SUBROUTINE CALL
!c      CALL NWFPEC(OSD,AO,L,IR,IL,N,MO,MJ0,MJ)
call nwfpec(aic,oaic,cv,ast1,osd,ao,aao,l,ir,il,n,mo,mj0,&
&lmax2,lcv1)
!     MO: MAICE AR-MODEL ORDER DETERMINED BY NWFPEC
!     RGT=R12*G' COMPUTATION
!c      CALL SBRUGT(MO,ID)
call sbrugt(mo,id,ast1,cv,rgt,mj1,ir,lmax2,lcv1)
mmmh=(mo+1)*id
ninew=mmmh
!     FRG' COMPUTATION
!     M1: NUMBER OF VARIABLE IN THE FUTURE SET
!     M2: NUMBER OF VARIABLES IN THE PAST SET =MMMH
m1=1
r11c=cv(1,1,1)
!-----------
!xx      DO 22 I=1,MJ1
!xx         DO 22 J=1,MJ1
!xx            FL(I,J)=CST0
!xx   22 CONTINUE
fl(1:mj1,1:mj1)=cst0
!------------
fl(1,1)=cst1/dsqrt(r11c)
fct=fl(1,1)
m2=(mo+1)*id
do 100 j=1,m2
frg(1,j)=fct*rgt(1,j)
100 continue
isw=0
ninew0=2
do 191 i=1,ir
!xx  191 NH(I)=I
nh(i)=i
191 continue
!xx      DO 634 I=1,MMMH
!xx  634 IH(I)=0
ih(1:mmmh)=0
500 m1=m1+1
m2=mmmh
m=m1+m2
mp1=m1+1
!     COVARIANCE MATRIX ARRANGEMENT
ninewi=ninew0
do 636 i=ninewi,mmmh
if(ih(i).eq.0) go to 637
!     TEST OF EXHAUSTION OF FUTURE VARIABLES
if(i.ge.mmmh) go to 1140
ninew0=ninew0+1
go to 636
637 nh(m1)=i
go to 638
636 continue
638 continue
!     THE AUGMENTED FRG' MATRIX (F+)(G+)G' IS OBTAINED BY THE
!     BORDERING TECHNIQUE:
!                           +-         -+
!                  (F+) = I F    0     I
!                           I             I
!                           I             I
!                           I SF'  SG  I
!                           +-         -+,
!                           +-         -+
!              (RVV+) = I RVV  SR I
!                           I             I
!                           I             I
!                           I SR'  SS  I
!                           +-        - +,
!              (F+)(RVV+)(F+)'=I,
!                           +-         -+
!                  (R+) = I    R      I
!                           I            I
!                           I            I
!                           I    ST'    I
!                           +-         -+,
!             AND
!                           +-         -+
!          (F+)(R+)G' = I   FRG'    I
!                           I             I
!                           I             I
!                           I  LRFRG  I
!                           +-         -+,
!        FRG' IS STORED AS FRG.
!     SR,SS,VTG=(ST)'G' ARRANGEMENT
ii=nh(m1)
idb=ii/id
ires=ii-idb*id
if(ires.ne.0) idb=idb+1
if(ires.eq.0) ires=id
do 640 j=1,m2
!xx  640 VTG(J)=RGT(II,J)
vtg(j)=rgt(ii,j)
640 continue
do 642 j=1,m1
jj=nh(j)
jdb=jj/id
jres=jj-jdb*id
if(jres.ne.0) jdb=jdb+1
if(jres.eq.0) jres=id
ijdb=idb-jdb
if(ijdb.lt.0) go to 646
sr(j)=cv(ijdb+1,ires,jres)
go to 642
646 ijdb1=-ijdb+1
sr(j)=cv(ijdb1,jres,ires)
642 continue
ss=sr(m1)
!     ZZ= F*SR COMPUTATION (F STORED AS FL IN COMMON AREA)
mm1=m1-1
!c      CALL BLMLVC(SR,ZZ,MM1)
call blmlvc(sr,zz,mm1,fl,mj1)
!     INNER PRODUCT OF FL*SR COMPUTATION
!     COMMON SUBROUTINE CALL
call innerp(zz,zz,fsrinp,mm1)
!     SG COMPUTATION
ssfr2=ss-fsrinp
!c      IF  (SSFR2) 731,731,732
if  (ssfr2 .lt. 0) go to 731
if  (ssfr2 .eq. 0) go to 731
if  (ssfr2 .gt. 0) go to 732
731 sg=cst0
go to 733
732 sg=cst1/dsqrt(ssfr2)
!     SF=-(F'*(F*SR))*SG COMPUTATION
!c  733 CALL AVMLVC(ZZ,SF,MM1)
733 call avmlvc(zz,sf,mm1,fl,mj1)
do 350 i=1,mm1
!xx  350 SF(I)=-SF(I)*SG
sf(i)=-sf(i)*sg
350 continue
!     F AUGMENTATION
do 360 j=1,mm1
!xx  360 FL(M1,J)=SF(J)
fl(m1,j)=sf(j)
360 continue
fl(m1,m1)=sg
!     SFRG=SF'*R12*G' COMPUTATION
!c      CALL VECMTX(SF,SFRG,NH,MM1,M2)
call vecmtx(sf,sfrg,nh,rgt,mm1,m2,mj1)
!      LRFRG=SFRG+SG*VTG (=THE LAST ROW OF (F+)(R+)G')
do 370 i=1,m2
!xx  370 FRG(M1,I)=SFRG(I)+SG*VTG(I)
frg(m1,i)=sfrg(i)+sg*vtg(i)
370 continue
!xx      DO 230 I=1,M2
do 231 i=1,m2
do 230 j=1,m1
u(i,j)=frg(j,i)
230 continue
231 continue
!     CANONICAL WEIGHTS FOR THE SET OF ORTHO-NORMALIZED FUTURE VARIABLES
!     ARE OBTAINED IN V(COMMON AREA).
!     CANONICAL CORRELATION COEFFICIENTS ARE RETURNED IN Z.
!     COMMON SUBROUTINE CALL
!c      CALL MSVD(Z,M2,M1)
nc = nc+1
call msvd(u,v,z(1,nc),m2,m1,mj1,mj1)
!     W=V'*FL : CANONICAL WEIGHTS COMPUTATION
!c      CALL MWTFL(M1)
call mwtfl(v,vv(1,1,nc),m1,fl,mj1)
do 260 j=1,m1
!c  260 Y(J)=Z(J)*Z(J)
!xx  260 Y(J,NC)=Z(J,NC)*Z(J,NC)
y(j,nc)=z(j,nc)*z(j,nc)
260 continue
!c      WRITE(6,6)
!c      WRITE(6,7) M1
!c      WRITE(6,8) M2
n1(nc)=m1
n2(nc)=m2
!c      WRITE(6,9) N
!c      WRITE(6,35)
!c      DO 8100 I=1,M1
!c 8100 WRITE(6,8200) I,(VV(I,J),J=1,M1)
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
!c      XX(J)=9999.0D-00
xx(j,nc)=9999.0d-00
go to 45
!c  145 XX(J)=-EN*DLOG(WL(J))
145 xx(j,nc)=-en*dlog(wl(j))
45 continue
!c      NDT=M1*M2
!c      ANDT=NDT
ndt(1,nc)=m1*m2
andt=ndt(1,nc)
!     DIC(=CHI-SQUARE-2.0*D.F.) COMPUTATION
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
!c   51 WRITE(6,50) J1,Z(J),Y(J),XX(J),NDT,X3(J)
ndt(j,nc)=(m1-j1)*(m2-j1)
andt=ndt(j,nc)
x3(j,nc)=xx(j,nc)-cst2*andt
51 continue
do 4300 j=2,m1
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
if(x3(m1,nc).gt.cst0) go to 110
6999 m1m=m1-1
if(m1m.le.0) go to 110
!     TRANSITION MATRIX (F) COMPUTATION
!c      AII=CST1/VV(M1,M1)
aii=cst1/vv(m1,m1,nc)
do 5100 i=1,m1m
iaw=iaw+1
!c      F(I)=-VV(M1,I)*AII
!c      VF(IAW)=F(I)
f(i,nc)=-vv(m1,i,nc)*aii
vf(iaw)=f(i,nc)
5100 continue
!c      WRITE(6,5200)
!     COMMON SUBROUTINE CALL
!c      CALL SUBVCP   (F,M1M)
m1=m1-1
isw=isw+1
do 1120 i=ninew0,mmmh,id
!xx 1120 IH(I)=1
ih(i)=1
1120 continue
if(isw.ge.id) go to 1100
110 if(ninew0.ge.ninew) go to 1100
ninew0=ninew0+1
go to 500
1140 m1=m1-1
!c 1100 WRITE(6,1130)
1100 continue
!     STRUCTURAL CHARACTERISTIC VECTOR PRINT AND PUNCH OUT
m1nh=m1
!c      WRITE(6,1131) (NH(I),I=1,M1)
!c      WRITE(7,1) IR,M1
!c      WRITE(7,1) (NH(I),I=1,M1)
!     F MATRIX (IN VECTOR FORM) PUNCH OUT
!c      WRITE(7,2) (VF(I),I=1,IAW)
!     INPUT MATRIX (G) COMPUTATION
!c 1011 CALL SUBBMA(AO,G,NH,M1,ID,MO,MJ0,MJ)
!xx 1011 CALL SUBBMA(AO,G,NH,M1,ID,MO,MJ1,MJ0)
call subbma(ao,g,nh,m1,id,mo,mj1,mj0)
!     INPUT MATRIX (G) PRINT AND PUNCH OUT
!c      WRITE(6,61)
!     COMMON SUBROUTINE CALL
!c      CALL SUBMPR(G,M1,ID,MJ1,MJ)
!c      IRP1=IR+1
!c      IF(IRP1.GT.M1) GO TO 1400
!c      DO 2010 I=IRP1,M1
!c 2010 WRITE(7,2) (G(I,J),J=1,ID)
!c 1400 CALL FLCLS3(NFL)
continue
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
!xx   50 FORMAT(/1H ,10X,I3,4X,F8.4,4X,F8.3,7X,F8.2,4X,I6,3X,F10.4)
!xx 4410 FORMAT(/1H ,'MINIMUM DIC(P) =',F12.2,1X,'ATTAINED AT P='
!xx     A,I5)
!xx 5200 FORMAT(/1H ,5X,'F(I)')
!xx 1130 FORMAT(/1H ,'STRUCTURAL CHARACTERISTIC VECTOR',
!xx     A' (H(I),I=1,P)')
!xx 1131 FORMAT(1H ,4X,10I12)
!xx   61 FORMAT(/1H ,'G-MATRIX')
!xx  259 FORMAT(1H ,6HINW(I),5X,10I5)
!xx 8200 FORMAT(/1H ,I5,10D12.5,/(1H ,5X,10D12.5))
!xx11111 FORMAT(/1H ,'INITIAL AUTO REGRESSIVE MODEL FITTING ',
!xx     A'BY THE MINIMUM AIC PROCEDURE.')
!xx11112 FORMAT(1H ,8X,'        THE VALUE OF CHI-SQUARE AND DIC (P) ',
!xx     A'CORRESPONDING TO CANONICAL R=1.000 ',
!xx     A'SHOULD BE IGNORED')
!xx11114 FORMAT(1H ,'PROGRAM 74.2.1. CANOCA')
!xx11113 FORMAT(/1H ,'IR=',I5)
end
!
!c      SUBROUTINE AVMLVC(Y,Z,MM)
subroutine avmlvc(y,z,mm,fl,mj1)
  use timsac_kinds, only: dp
  implicit none
!     Z=FL'*Y
!     FL: LOWER TRIANGLE
!     Y: VECTOR
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION Y(MM),Z(MM)
!c      DIMENSION FL(46,46)
!c      COMMON /COM25/FL
!xx      DIMENSION FL(MJ1,MJ1)
integer mm ,mj1
real(dp) y(mm), z(mm), fl(mj1,mj1)
! local
integer i, j
real(dp) cst0, sum
cst0=0.0d-00
do 10 i=1,mm
sum=cst0
do 11 j=i,mm
!xx   11 SUM=SUM+FL(J,I)*Y(J)
sum=sum+fl(j,i)*y(j)
11 continue
z(i)=sum
10 continue
return
end
!
!c      SUBROUTINE BLMLVC(Y,Z,MM)
subroutine blmlvc(y,z,mm,fl,mj1)
  use timsac_kinds, only: dp
  implicit none
!     Z=FL*Y
!     FL: LOWER TRIANGLE
!     Y: VECTOR
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION Y(MM),Z(MM)
!c      DIMENSION FL(46,46)
!c      COMMON /COM25/FL
!xx      DIMENSION FL(MJ1,MJ1)
integer mm, mj1
real(dp) y(mm), z(mm), fl(mj1,mj1)
! local
integer i, j
real(dp) cst0, sum
cst0=0.0d-00
do 10 i=1,mm
sum=cst0
do 11 j=1,i
!xx   11 SUM=SUM+FL(I,J)*Y(J)
sum=sum+fl(i,j)*y(j)
11 continue
z(i)=sum
10 continue
return
end
!
!c      SUBROUTINE MWTFL(V,VV,MM)
subroutine mwtfl(v,vv,mm,fl,mj1)
  use timsac_kinds, only: dp
  implicit none
!     VV=V'*FL
!xx      IMPLICIT REAL *8(A-H,O-Z)
!c      COMMON /COM20/V
!c      COMMON /COM21/VV
!c      COMMON /COM25/FL
!c      DIMENSION V(46,46),VV(46,46),FL(46,46)
!xx      DIMENSION V(MJ1,MJ1),VV(MJ1,MJ1),FL(MJ1,MJ1)
integer mm, mj1
real(dp) v(mj1,mj1), vv(mj1,mj1), fl(mj1,mj1)
! local
integer i, j, k
real(dp) cst0, sum
cst0=0.0d-00
do 10 i=1,mm
do 11 j=1,mm
sum=cst0
do 12 k=1,mm
!xx   12 SUM=SUM+V(K,I)*FL(K,J)
sum=sum+v(k,i)*fl(k,j)
12 continue
vv(i,j)=sum
11 continue
10 continue
return
end
!
!c      SUBROUTINE NLTIV(R,RIN,DET,K,MJ)
subroutine nltiv(r,rin,det,k)
  use timsac_kinds, only: dp
  implicit none
!     INVERSE OF R IS FACTORI
!     INVERSE OF R IS FACTORED INTO THE FORM L'*L IS
!     INVERSE OF R IS FACTORED INTO THE FORM L'*L. L,LOWER
!     INVERSE OF R IS FACTORED INTO THE FORM L'*L. L, LOWER TRIANGULAR,
!     IS RETURNED IN R. OTHER OUTPUTS ARE
!     RIN: INVERSE OF DIAGONAL L(I,J)
!     DET: DETERMINANT OF R(I,J)
!     MJ IS THE ABSOLUTE DIMENSION OF R IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION R(MJ,MJ),RIN(K)
!xx      DIMENSION R(K,K),RIN(K)
integer k
real(dp) r(k,k), rin(k), det
! local
integer i, km1, l, l1, m
real(dp) cst0, cst1, abc0, rr, rpivot, ril
cst0=0.0d-00
cst1=1.0d-00
km1=k-1
abc0=cst0
det=cst1
do 10 l=1,k
rr=r(l,l)
det=det*rr
rpivot=cst1/dsqrt(rr)
r(l,l)=rpivot
rin(l)=cst1/rpivot
do 12 i=1,k
if(i.eq.l) go to 12
r(l,i)=rpivot*r(l,i)
12 continue
if(l.eq.k) go to 11
l1=l+1
do 13 i=l1,k
ril=-rpivot*r(i,l)
r(i,l)=ril*rpivot
do 14 m=1,k
if(m.eq.l) go to 14
r(i,m)=r(i,m)+ril*r(l,m)
14 continue
13 continue
10 continue
11 return
end
!
!c      SUBROUTINE NWFPEC(OSD,AO,L,IR,IL,N,IFPEC,MJ0,MJ)
subroutine nwfpec(aaic,oaic,cv,ast1,osd,ao,aao,l,ir,il,n,ifpec,&
&mj0,lmax2,lcv1)
  use timsac_kinds, only: dp
  implicit none
!     AR-FITTING
!     AUTOREGRESSIVE MODEL FITTING BY THE MINIMUM AIC PROCEDURE.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION CV(25,7,7)
!c      DIMENSION A1(13,7,7),B1(13,7,7)
!c      DIMENSION SD(7,7),SE(7,7),SF(7,7)
!c      DIMENSION XSD(7,7),XSF(7,7),D(7,7),E(7,7),Z1(7,7)
!c      DIMENSION SFL(7,7),RIN(7),B(7,7)
!c      DIMENSION AST1(91,7,7)
!c      DIMENSION OSD(7,7),AO(13,7,7)
!xx      DIMENSION CV(LCV1,IR+IL,IR+IL)
!xx      DIMENSION A1(L,IR+IL,IR+IL),B1(L,IR+IL,IR+IL)
!xx      DIMENSION SD(IR+IL,IR+IL),SE(IR+IL,IR+IL),SF(IR+IL,IR+IL)
!xx      DIMENSION XSD(IR+IL,IR+IL),XSF(IR+IL,IR+IL)
!xx      DIMENSION D(IR+IL,IR+IL),E(IR+IL,IR+IL),Z1(IR+IL,IR+IL)
!xx      DIMENSION SFL(IR+IL,IR+IL),RIN(IR+IL),B(IR+IL,IR+IL)
!xx      DIMENSION AST1(LMAX2,IR+IL,IR+IL)
!xx      DIMENSION OSD(IR,IR),AO(MJ0,IR,IR+IL),AAO(MJ0,IR,IR+IL)
!xx      DIMENSION AAIC(0:MJ0-1)
integer l, ir, il, n, ifpec, mj0, lmax2, lcv1
real(dp) aaic(0:mj0-1), oaic, cv(lcv1,ir+il,ir+il),&
&ast1(lmax2,ir+il,ir+il), osd(ir,ir),&
&ao(mj0,ir,ir+il), aao(mj0,ir,ir+il)
! local
integer i, ii, imi, inx, ip, jj, l1, m, ms
real(dp) a1(l,ir+il,ir+il), b1(l,ir+il,ir+il),&
&sd(ir+il,ir+il), se(ir+il,ir+il),&
&sf(ir+il,ir+il), xsd(ir+il,ir+il),&
&xsf(ir+il,ir+il), d(ir+il,ir+il), e(ir+il,ir+il),&
&z1(ir+il,ir+il), sfl(ir+il,ir+il), rin(ir+il),&
&b(ir+il,ir+il), aic, sfdt, sddet, sfdet, cst0
!c      COMMON /COM9/AST1
!c      COMMON /COM10/CV
!
ip=ir+il
l1=l+1
cst0=0.0d-00
!     INITIAL CONDITION AND COVARIANCE PRINT OUT
!     INITIAL SD, SF, SE COMPUTATION
!xx      DO 330 II=1,IP
do 331 ii=1,ip
do 330 jj=1,ip
sd(ii,jj)=cv(1,ii,jj)
sf(ii,jj)=sd(ii,jj)
sfl(ii,jj)=sf(ii,jj)
se(ii,jj)=cv(2,ii,jj)
xsd(ii,jj)=sd(ii,jj)
!xx  330 XSF(II,JJ)=SF(II,JJ)
xsf(ii,jj)=sf(ii,jj)
330 continue
331 continue
!     0-TH STEP COMPUTATION
ifpec=0
ms=0
!     AIC COMPUTATION
!c      CALL SAIC (SD,N,IP,MS,AIC,MJ)
call saic (sd,n,ip,ms,aic)
oaic=aic
!     AIC PRINT OUT
!c      WRITE(6,600)
!c      WRITE(6,264) MS,AIC
aaic(0)=aic
!c      CALL NLTIV(SFL,RIN,SFDT,IP,MJ)
call nltiv(sfl,rin,sfdt,ip)
!     SFL IS FACTORED INTO L'*L AND L IS RETURNED IN SFL
inx=1
!xx      DO 8 II=1,IP
do 9 ii=1,ip
do 8 jj=1,ii
ast1(1,jj,ii)=cst0
!xx    8 AST1(1,II,JJ)=SFL(II,JJ)
ast1(1,ii,jj)=sfl(ii,jj)
8 continue
9 continue
!     ITERATION M=1 TO L
do 400 m=1,l
!     INVERSE OF SD, SF COMPUTATION
!     COMMON SUBROUTINE CALL
!c      CALL INVDET(XSD,SDDET,IP,MJ)
call invdet(xsd,sddet,ip,ip)
!     COMMON SUBROUTINE CALL
!c      CALL INVDET(XSF,SFDET,IP,MJ)
call invdet(xsf,sfdet,ip,ip)
!     D, E, SD, SF COMPUTATION
!     COMMON SUBROUTINE CALL
!c      CALL MULPLY(SE,XSF,D,IP,IP,IP,MJ,MJ,MJ)
call mulply(se,xsf,d,ip,ip,ip)
!     COMMON SUBROUTINE CALL
!c      CALL TRAMDL(SE,XSD,E,IP,IP,IP,MJ,MJ,MJ)
call tramdl(se,xsd,e,ip,ip,ip)
!     COMMON SUBROUTINE CALL
!c      CALL TRAMDR(D,SE,Z1,IP,IP,IP,MJ,MJ,MJ)
call tramdr(d,se,z1,ip,ip,ip)
!     COMMON SUBROUTINE CALL
!c      CALL SUBTAL(SD,Z1,IP,IP,MJ,MJ)
call subtal(sd,z1,ip,ip)
!     COMMON SUBROUTINE CALL
!c      CALL MULPLY(E,SE,Z1,IP,IP,IP,MJ,MJ,MJ)
call mulply(e,se,z1,ip,ip,ip)
!     COMMON SUBROUTINE CALL
!c      CALL SUBTAL(SF,Z1,IP,IP,MJ,MJ)
call subtal(sf,z1,ip,ip)
ms=m
!xx      DO 410 II=1,IP
do 411 ii=1,ip
do 410 jj=1,ip
xsd(ii,jj)=sd(ii,jj)
xsf(ii,jj)=sf(ii,jj)
sfl(ii,jj)=sf(ii,jj)
410 continue
411 continue
!     AIC COMPUTATION
!c     CALL SAIC (SD,N,IP,MS,AIC,MJ)
call saic (sd,n,ip,ms,aic)
!     AIC PRINT OUT
!c     WRITE(6,264) MS,AIC
aaic(ms) = aic
!c      CALL NLTIV(SFL,RIN,SFDT,IP,MJ)
call nltiv(sfl,rin,sfdt,ip)
!     SFL IS FACTORED INTO L'*L AND L IS RETURNED IN SFL
!     FORWARD AND BACKWARD PREDICTOR COMPUTATION
!c      CALL COEFAB(A1,B1,D,E,MS,IP,MJ0,MJ)
call coefab(a1,b1,d,e,ms,l,ip)
do 14 i=1,m
imi=m+1-i
!xx      DO 15 II=1,IP
do 155 ii=1,ip
do 15 jj=1,ip
!xx   15 B(II,JJ)=B1(IMI,II,JJ)
b(ii,jj)=b1(imi,ii,jj)
15 continue
155 continue
!c      CALL BLMULP(SFL,B,Z1,IP,IP,MJ,MJ)
call blmulp(sfl,b,z1,ip,ip)
inx=inx+1
!xx      DO 16 II=1,IP
do 166 ii=1,ip
do 16 jj=1,ip
!xx   16 AST1(INX,II,JJ)=-Z1(II,JJ)
ast1(inx,ii,jj)=-z1(ii,jj)
16 continue
166 continue
14 continue
inx=inx+1
!xx      DO 17 II=1,IP
do 177 ii=1,ip
do 17 jj=1,ii
ast1(inx,jj,ii)=cst0
!xx   17 AST1(INX,II,JJ)=SFL(II,JJ)
ast1(inx,ii,jj)=sfl(ii,jj)
17 continue
177 continue
!     MINIMUM AIC SEARCH
if(oaic.le.aic) go to 440
oaic=aic
ifpec=m
!xx      DO 560 II=1,IR
do 563 ii=1,ir
do 560 jj=1,ir
!xx  560 OSD(II,JJ)=SD(II,JJ)
osd(ii,jj)=sd(ii,jj)
560 continue
563 continue
do 561 i=1,m
!xx      DO 562 II=1,IR
do 564 ii=1,ir
do 562 jj=1,ip
aao(i,ii,jj)=-a1(i,ii,jj)
!xx  562 AO(I,II,JJ)=-A1(I,II,JJ)
ao(i,ii,jj)=-a1(i,ii,jj)
562 continue
564 continue
561 continue
440 if(m.eq.l) go to 400
!     SE COMPUTATION
!     COMMON SUBROUTINE CALL
!c      CALL NEWSE(A1,SE,MS,IP,MJ0,MJ)
call newse(a1,cv,se,ms,l,ip,lcv1)
400 continue
!     MIN.AIC PRINT OUT
!c      WRITE(6,607) OAIC,IFPEC
!     OSD, AO PRINT OUT
!c      WRITE(6,608)
!     COMMON SUBROUTINE CALL
!c      CALL SUBMPR(OSD,IR,IR,MJ,MJ)
if(ifpec.le.0) go to 699
!c      WRITE(6,1601)
!c      WRITE(6,1602)
!c      WRITE(6,609)
!     COMMON SUBROUTINE CALL
!c      CALL PRMAT3(AO,IFPEC,IR,IP,0,MJ0,MJ,MJ)
do 1611 i=1,m
do 1612 ii=1,ir
!xx      DO 1612 JJ=1,IP
do 1613 jj=1,ip
!xx 1612 AO(I,II,JJ)=-AO(I,II,JJ)
ao(i,ii,jj)=-ao(i,ii,jj)
1613 continue
1612 continue
1611 continue
699 return
!xx    1 FORMAT(10I5)
!xx    2 FORMAT(4D20.10)
!xx   42 FORMAT(//1H ,17HCOVARIANCE MATRIX)
!xx  264 FORMAT(1H ,I5,2X,D16.7)
!xx  600 FORMAT(///1H ,4X,1HI,15X,3HAIC)
!xx  607 FORMAT(/1H ,'MINIMUM AIC ',D12.5,2X,'ATTAINED AT M=',I5)
!xx  608 FORMAT(//1H ,10X,10HOSD(II,JJ),': INNOVATION VARIANCE')
!xx  609 FORMAT(//1H ,10X,'A(I): AUTOREGRESSIVE COEFFICIENTS')
!xx 1601 FORMAT(/1H ,'AR-MODEL:')
!xx 1602 FORMAT(1H ,'Y(N)+A(1)Y(N-1)+...+A(L)Y(N-L)=X(N)')
end
!
!c      SUBROUTINE SAIC (SD,N,K,MS,AIC,MJ)
subroutine saic (sd,n,k,ms,aic)
  use timsac_kinds, only: dp
  implicit none
!     AIC COMPUTATION.
!     SD: COVARIANCE MATRIX OF INNOVATION
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION SD(MJ,MJ)
!c      DIMENSION SD1(7,7)
!xx      DIMENSION SD(K,K)
!xx      DIMENSION SD1(K,K)
integer n, k, ms
real(dp) sd(k,k), aic
! local
integer i, j
real(dp) sd1(k,k), an, sdrm, arm2
an=n
!xx      DO 9 I=1,K
do 10 i=1,k
do 9 j=1,k
!xx    9 SD1(I,J)=SD(I,J)
sd1(i,j)=sd(i,j)
9 continue
10 continue
!     COMMON SUBROUTINE CALL
!c      CALL SUBDET(SD1,SDRM,K,MJ)
call subdetc(sd1,sdrm,k)
arm2=2*ms*k*k
aic=an*dlog(sdrm)+arm2
return
end
!
!c      SUBROUTINE SBRUGT(MO,ID)
subroutine sbrugt(mo,id,ast1,cv,rgt,mj1,mj,lmax2,lcv1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES MATRIX R12*G'.
!     INPUTS REQUIRED:
!     ID: DIMENSION OF Y(I) =IR=IP
!     MO: MAICE ORDER OF AR-MODEL
!     CV(I,II,JJ): COVARIANCE MATRIX
!     AST1(I,II,JJ): MATRIX G FOR ORTHO-NORMALIZATION OF THE PRESENT AND
!     PAST VARIABLES. G IS A LOWER TRIANGULAR BLOCK MATRIX WITH THE
!     STRUCTURE
!                  +-                            -+
!             G = I  L0	     0        0   ....  I
!                  I -L1B11  L1       0  .....  I
!                  I -L2B22 -L2B21 L2 ....  I
!                  I   .        .        ..      	...  I
!                  I   .        .        . .	..   I
!                  I   .        .        .  .	.    I
!                  I   .        .	   .   .       I
!                  I   .        .	   .    .      I
!                  I   .        .	   .     .     I
!                  I   .        .	   .      .    I
!                  +-			   -+
!     WHERE BMI'S ARE THE COEFFICIENTS OF THE M-TH ORDER BACKWARD AUTO-
!     REGRESSION,
!     Y(I-M)-BM1Y(I-M+1)-..-BMMY(I)=ZM(I-M)
!     AND LM(LM)'=INVERSE OF SD(M),WHERE SD(M) IS THE COVARIANCE MATRIX
!     RESIDUAL ZM(I).
!     OUTPUT:
!     RGT(I,J): MATRIX OF R12*G'
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION CV(25,7,7),AST1(91,7,7)
!c      DIMENSION RGT(91,91)
!c      DIMENSION X(7,7),Y(7,7)
!xx      DIMENSION CV(LCV1,ID,ID),AST1(LMAX2,MJ,MJ)
!xx      DIMENSION RGT(MJ1,MJ1)
!xx      DIMENSION X(ID,ID),Y(ID,ID)
integer mo, id, mj1, mj, lmax2, lcv1
real(dp) ast1(lmax2,mj,mj), cv(lcv1,id,id), rgt(mj1,mj1)
! local
integer i, ib, ii, im1, in, inc, inx, irg, j, jb, jj, jn, jnc, k,&
&kk, mj7, mp1
real(dp) x(id,id), y(id,id), cst0, sum
!      COMMON /COM9/AST1
!c      COMMON /COM10/CV
!c      COMMON /COM11/RGT
cst0=0.0d-00
mp1=mo+1
!c      IRG=91
irg=mj1
!xx      DO 9 I=1,IRG
!xx      DO 9 J=1,IRG
!xx      RGT(I,J)=CST0
!xx    9 CONTINUE
rgt(1:irg,1:irg)=cst0
inc=0
do 10 i=1,mp1
im1=i-1
inx=0
jnc=0
do 11 j=1,mp1
do 12 k=1,j
ib=im1+k
jb=inx+k
do 18 ii=1,id
!xx      DO 18 JJ=1,ID
do 19 jj=1,id
x(ii,jj)=cv(ib,ii,jj)
y(ii,jj)=ast1(jb,ii,jj)
19 continue
18 continue
mj7=7
do 15 ii=1,id
do 16 jj=1,id
sum=cst0
do 17 kk=1,id
!xx   17 SUM=SUM+X(II,KK)*Y(JJ,KK)
sum=sum+x(ii,kk)*y(jj,kk)
17 continue
in=inc+ii
jn=jnc+jj
rgt(in,jn)=rgt(in,jn)+sum
16 continue
15 continue
12 continue
jnc=jnc+id
inx=inx+j
11 continue
inc=inc+id
10 continue
return
end
!
!c      SUBROUTINE SUBBMA(AO,B,NH,M1,ID,IQ,MJ0,MJ)
subroutine subbma(ao,b,nh,m1,id,iq,mj1,mj0)
  use timsac_kinds, only: dp
  implicit none
!     AR-FITTING
!     B-MATRIX COMPUTATION
!     (USE AR-COEFFICIENTS)
!     M1 IS LESS THAN 47.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION AO(MJ0,MJ,MJ)
!c      DIMENSION X(13,7,7)
!c      DIMENSION XX(7,7),AA(7,7),C(7,7)
!c      DIMENSION W(100,7),B(46,7)
!c      DIMENSION NH(91)
!xx      DIMENSION AO(MJ0,ID,ID)
!xx      DIMENSION X(IQ,ID,ID)
!xx      DIMENSION XX(ID,ID),AA(ID,ID),C(ID,ID)
!xx      DIMENSION W(100,ID),B(MJ1,ID)
!xx      DIMENSION NH(M1)

integer m1, id, iq, mj1, mj0, nh(m1)
real(dp) ao(mj0,id,id), b(mj1,id)
! local
integer i, ic, ii, iim1, ij, j, jj, mq1, nhc
real(dp) x(iq,id,id), xx(id,id), aa(id,id), c(id,id),&
&w(100,id), cst0, cst1
cst0=0.0d-00
cst1=1.0d-00
w(1:id,1:id)=cst0
do 10 i=1,id
!xx      DO 11 J=1,ID
!xx   11 W(I,J)=CST0
!xx   10 W(I,I)=CST1
w(i,i)=cst1
10 continue
if(iq.le.1) go to 110
mq1=iq-1
do 200 ii=1,mq1
!xx      DO 210 I=1,ID
do 211 i=1,id
do 210 j=1,id
!xx  210 X(II,I,J)=AO(II,I,J)
x(ii,i,j)=ao(ii,i,j)
210 continue
211 continue
if(ii.le.1) go to 201
iim1=ii-1
do 220 jj=1,iim1
!xx      DO 230 I=1,ID
do 231 i=1,id
do 230 j=1,id
!xx  230 AA(I,J)=AO(JJ,I,J)
aa(i,j)=ao(jj,i,j)
230 continue
231 continue
ij=ii-jj
!xx      DO 240 I=1,ID
do 241 i=1,id
do 240 j=1,id
!xx  240 XX(I,J)=X(IJ,I,J)
xx(i,j)=x(ij,i,j)
240 continue
241 continue
!c      CALL MULPLY(AA,XX,C,ID,ID,ID,MJ,MJ,MJ)
call mulply(aa,xx,c,id,id,id)
!xx      DO 250 I=1,ID
do 251 i=1,id
do 250 j=1,id
!xx  250 X(II,I,J)=X(II,I,J)+C(I,J)
x(ii,i,j)=x(ii,i,j)+c(i,j)
250 continue
251 continue
220 continue
!xx  201 DO 260 I=1,ID
201 do 261 i=1,id
do 260 j=1,id
ic=id+(ii-1)*id+i
!xx  260 W(IC,J)=X(II,I,J)
w(ic,j)=x(ii,i,j)
260 continue
261 continue
200 continue
110 do 300 i=1,m1
nhc=nh(i)
do 310 j=1,id
!xx  310 B(I,J)=W(NHC,J)
b(i,j)=w(nhc,j)
310 continue
300 continue
!xx  100 RETURN
return
end
!
!c      SUBROUTINE VECMTX(X,Z,NH,MM,NN)
subroutine vecmtx(x,z,nh,rgt,mm,nn,mj1)
  use timsac_kinds, only: dp
  implicit none
!     Z=X*Y (X,Z: VECTORS, Y: SUBMATRIX)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION X(MM),Z(NN)
!c      DIMENSION RGT(91,91)
!c      DIMENSION NH(91)
!xx      DIMENSION RGT(MJ1,MJ1)
!xx      DIMENSION NH(MM)
integer mm, nn, mj1, nh(mm)
real(dp) x(mm), z(nn), rgt(mj1,mj1)
! local
integer i, j, jj
real(dp) cst0, sum
!c      COMMON /COM11/RGT
cst0=0.0d-00
do 10 i=1,nn
sum=cst0
do 11 j=1,mm
jj=nh(j)
!xx   11 SUM=SUM+X(J)*RGT(JJ,I)
!xx   10 Z(I)=SUM
sum=sum+x(j)*rgt(jj,i)
11 continue
z(i)=sum
10 continue
return
end
!
!c      SUBROUTINE BLMULP(X,Y,Z,MM,NN,MJ1,MJ2)
subroutine blmulp(x,y,z,mm,nn)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE
!     Z=X*Y (X: LOWER TRIANGLE)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION X(MJ1,MJ1),Y(MJ1,MJ2),Z(MJ1,MJ2)
!xx      DIMENSION X(MM,MM),Y(MM,NN),Z(MM,NN)
integer mm, nn
real(dp) x(mm,mm), y(mm,nn), z(mm,nn)
! local
integer i, j, k
real(dp) cst0, sum
cst0=0.0d-00
do 10 i=1,mm
do 11 j=1,nn
sum=cst0
do 12 k=1,i
!xx   12 SUM=SUM+X(I,K)*Y(K,J)
sum=sum+x(i,k)*y(k,j)
12 continue
z(i,j)=sum
11 continue
10 continue
return
end
!
!
!c      SUBROUTINE SUBDET(X,XDETMI,MM,MJ)
subroutine subdetc(x,xdetmi,mm)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE
!     THIS SUBROUTINE COMPUTES THE DETERMINANT OF UPPER LEFT MM X MM
!     OF X.  FOR GENERAL USE STATEMENTS 20-21 SHOULD BE RESTORED.
!     X: ORIGINAL MATRIX
!     XDETMI: DETERMINANT OF UPPER LEFT MM X MM OF X
!     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(X)
!c      DIMENSION X(MJ,MJ)
!xx      DIMENSION X(MM,MM)
integer mm
real(dp) x(mm,mm), xdetmi
! local
integer i, i1, j, k, mm1
real(dp) cst0, cst1, xc, xxc
cst0=0.0d-00
cst1=1.0d-00
xdetmi=cst1
if(mm.eq.1) go to 18
mm1=mm-1
do 10 i=1,mm1
!   20 IF(X(I,I).NE.CST0) GO TO 11
!     DO 12 J=I,MM
!     IF(X(I,J).EQ.CST0) GO TO 12
!     JJ=J
!     GO TO 13
!  12 CONTINUE
!     XDETMI=CST0
!     GO TO 17
!  13 DO 14 K=I,MM
!     XXC=X(K,JJ)
!     X(K,JJ)=X(K,I)
!  14 X(K,I)=XXC
!  21 XDETMI=-XDETMI
!xx   11 XDETMI=XDETMI*X(I,I)
xdetmi=xdetmi*x(i,i)
xc=cst1/x(i,i)
i1=i+1
do 15 j=i1,mm
xxc=x(j,i)*xc
do 16 k=i1,mm
!xx   16 X(J,K)=X(J,K)-X(I,K)*XXC
x(j,k)=x(j,k)-x(i,k)*xxc
16 continue
15 continue
10 continue
18 xdetmi=xdetmi*x(mm,mm)
!xx   17 RETURN
return
end
