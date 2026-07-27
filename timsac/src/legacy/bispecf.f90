! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine bispecf(n,mh,cc0,c0,p1,p2,q,a,br,bi,rat)
  use timsac_kinds, only: dp
  implicit none
!
!c	PROGRAM BISPEC
!     PROGRAM 74.6.2.
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(2
!         BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
!         NO.6 MARCH 1976, THE INSTITUTE OF STATISTICAL MATHEMATICS
!-----------------------------------------------------------------------
!     THIS PROGRAM COMPUTES BISPECTRUM USING THE DIRECT FOURIER TRANSFOR
!     OF SAMPLE THIRD ORDER MOMENTS.
!     THIS PROGRAM REQUIRES THE FOLLOWING INPUTS;
!     OUTPUTS OF THE PROGRAM THIRMO:
!      N; DATA LENGTH,
!      MH; MAXIMUM LAG,
!      CC(I); AUTOCOVARIANCES,
!      C(I,J); THIRD ORDER MOMENTS.
!
!c      !DEC$ ATTRIBUTES DLLEXPORT :: BISPECF
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION C(51,51),C1(102),S1(102)
!c      DIMENSION CL(51,51),SL(51,51)
!c      DIMENSION CA1(51,51),CA2(51,51)
!c      DIMENSION P(51)
!c      DIMENSION CC(51)
!xx      DIMENSION C(MH+1,MH+1),C1((MH+1)*2),S1((MH+1)*2)
!xx      DIMENSION CL(MH+1,MH+1),SL(MH+1,MH+1)
!xx      DIMENSION CA1(MH+1,MH+1),CA2(MH+1,MH+1)
!xx      DIMENSION P1(MH+1),P2(MH+1),Q(MH+1)
!xx      DIMENSION CC(MH+1)
!xx      DIMENSION CC0(MH+1),C0(MH+1,MH+1)
!xx      DIMENSION A(MH+1,MH+1),BR(MH+1,MH+1),BI(MH+1,MH+1)
integer n, mh
real(dp) cc0(mh+1), c0(mh+1,mh+1), p1(mh+1), p2(mh+1),&
&q(mh+1), a(mh+1,mh+1), br(mh+1,mh+1),&
&bi(mh+1,mh+1), rat
! local
integer i, i2h, id, id1, ir, ir1, irl, is, is1, j, l1, mc1, mc2,&
&mc3, md, mdr, mds, mh1, mms, mms1, mr, mrd, mrs, ms, msd,&
&msr, ns, ns1
real(dp) c(mh+1,mh+1), c1((mh+1)*2), s1((mh+1)*2),&
&cl(mh+1,mh+1), sl(mh+1,mh+1), ca1(mh+1,mh+1),&
&ca2(mh+1,mh+1), cc(mh+1), cst0, cst1, cst2,&
&cst6, cst2i, cst6i, h, pi, pp, ai, t, c00, sl0,&
&tc1, tc2, tc3, ts1, ts3, cn0, sn0
!
!     INPUT / OUTPUT DATA FILE OPEN
!c	CALL SETWND
!c	CALL FLOPN2(NFL)
!c	IF (NFL.EQ.0) GO TO 999
!
mh1 = mh+1
a(1:mh1,1:mh1) = 0.0d-00
br(1:mh1,1:mh1) = 0.0d-00
bi(1:mh1,1:mh1) = 0.0d-00
cl(1:mh1,1:mh1) = 0.0d-00
sl(1:mh1,1:mh1) = 0.0d-00
!
!c      MJ1=51
cst0 =0.0d-00
cst1 =1.0d-00
cst2 =2.0d-00
cst6=6.0d-00
cst2i=1.0d-00/cst2
cst6i=1.0d-00/cst6
!     INITIAL LOADING
!     AUTOCOVARIANCE AND THIRD ORDER MOMENT INPUT
!c	READ(5,1) N,MH
l1=mh+1
!c	READ(5,2) (CC(I),I=1,L1)
h=mh
!xx      DO 10 I=1,L1
do 11 i=1,l1
!c   10 READ(5,2) (C(I,J),J=1,I)
cc(i)=cc0(i)
do 10 j=1,i
c(i,j)=c0(i,j)
10 continue
11 continue
!c      WRITE(6,130)
!c      WRITE(6,129) N,MH
!     AUTOCOVARIANCES PRINT OUT
!c      WRITE(6,163)
!c      DO 164 I=1,L1
!c      IM1=I-1
!c  164 WRITE(6,165) IM1,CC(I)
!     POWER SPECTRUM COMPUTATION
!c      CALL SAUSP1(CC,P,N,L1,L1)
call sausp1(cc,p1,p2,q,n,l1,l1)
!     THIRD ORDER MOMENTS PRINT OUT
!c      WRITE(6,180)
!c      CALL MOMTPR (C,MJ1,L1)
!     END CORRECTIONS
c(1,1)=c(1,1)*cst6i
do 112 i=2,l1
c(i,1) =c(i,1)*cst2i
c(i,i)=c(i,i)*cst2i
!xx  112 C(L1,I)=C(L1,I)*CST2I
c(l1,i)=c(l1,i)*cst2i
112 continue
c(l1,1) =c(l1,1)*cst2i
!     FOURIER TRANSFORMATION
pi=3.141592653
pp=pi/h
i2h=mh+mh
do 115 i=1,i2h
ai=i
t=pp*ai
c1(i)=dcos(t)
!xx  115 S1(I)=DSIN(T)
s1(i)=dsin(t)
115 continue
!     CL(0,0) COMPUTATION
c00 =cst6*c(1,1)
sl0 =cst0
is=0
id=0
ir=0
tc1 =cst0
do 600 mms1=1,l1
tc2=cst0
do 610 ns1=1,mms1
!xx  610 TC2=TC2+C(MMS1,NS1)
tc2=tc2+c(mms1,ns1)
610 continue
tc1=tc1+tc2
600 continue
cl(1,1)=cst6*tc1
!     CL(IS,0) COMPUTATION
do 120 is=1,mh
id=is
tc3=cst0
ms=0
do 122 mms=1,mh
ms=ms+is
if(ms.le.i2h) go to 461
ms=ms-i2h
461 mc1=ms
mc2=ms
mc3=0
tc2 =c(mms+1,1)*(c1(mc1)+c1(mc2)+cst1)
tc1 =cst0
do 123 ns=1,mms
mc2=mc2-is
mc3=mc3+is
if(mc2.gt.0) go to 561
mc2=mc2+i2h
561 if(mc3.le.i2h) go to 562
mc3=mc3-i2h
562 tc1=tc1+c(mms+1,ns+1)*(c1(mc1)+c1(mc2)+c1(mc3))
123 continue
tc3=tc3+tc2+tc1
122 continue
is1=is+1
cl(is1,1)=c00+tc3+tc3
cl(1,is1)=cl(is1,1)
sl(is1,1) =cst0
sl(1,is1)=sl(is1,1)
120 continue
!     CL(ID,IR),SL(ID,IR) COMPUTATION
do 20 is=2,mh
irl=is/2
do 21 ir=1,irl
tc3=cst0
ts3 =cst0
id=is-ir
md=0
ms=0
mr=0
do 22 mms=1,mh
md=md+id
ms=ms-is
mr=mr+ir
if(md.le.i2h) go to 401
md=md-i2h
401 if(ms.gt.0) go to 412
ms=ms+i2h
412 if(mr.le.i2h) go to 413
mr=mr-i2h
413 mdr=md
mds=mdr
msr=ms
msd=msr
mrd=mr
mrs=mrd
cn0=c(mms+1,1)*(c1(mdr)+c1(mds)+c1(msr)+c1(msd)+c1(mrd)+c1(mrs))
sn0=c(mms+1,1)*(s1(mdr)+s1(mds)+s1(msr)+s1(msd)+s1(mrd)+s1(mrs))
tc1 =cst0
ts1=cst0
do 23 ns=1,mms
msd=msd+id
mrd=mrd+id
mdr=mdr+ir
msr=msr+ir
mds=mds-is
mrs=mrs-is
if(msd.le.i2h) go to 501
msd=msd-i2h
501 if(mrd.le.i2h) go to 511
mrd=mrd-i2h
511 if(mdr.le.i2h) go to 521
mdr=mdr-i2h
521 if(msr.le.i2h) go to 531
msr=msr-i2h
531 if(mds.gt.0) go to 542
mds=mds+i2h
542 if(mrs.gt.0) go to 552
mrs=mrs+i2h
552 tc1=tc1+c(mms+1,ns+1)*(c1(msd)+c1(mrd)+c1(mdr)+c1(msr)+c1(mds)+c1(&
&mrs))
ts1=ts1+c(mms+1,ns+1)*(s1(msd)+s1(mrd)+s1(mdr)+s1(msr)+s1(mds)+s1(&
&mrs))
23 continue
tc3=tc3+cn0+tc1
ts3=ts3+sn0+ts1
22 continue
ir1=ir+1
id1=id+1
cl(id1,ir1)=c00+tc3
sl(id1,ir1)=ts3
cl(ir1,id1)=cl(id1,ir1)
sl(ir1,id1)=sl(id1,ir1)
21 continue
20 continue
!     SMOOTHING OF THE RAW ESTIMATES
!     ROWWISE SMOOTHING
call subca(cl,ca1,mh,0)
call subca(sl,ca2,mh,1)
!     COLUMNWISE SMOOTHING
call subcb(ca1,cl,mh)
call subcb(ca2,sl,mh)
!     OBLIQUE SMOOTHING
!c      WRITE(6,160)
!c      CALL SUBCD(CL,CA1,MH)
!c      WRITE(6,161)
!c      CALL SUBCD(SL,CA2,MH)
!     Q1(I,J) COMPUTATION
!c      CALL SUBQ1(CA1,CA2,P,N,MH)
call subcd(cl,ca1,mh,br)
call subcd(sl,ca2,mh,bi)
call subq1(ca1,ca2,p1,n,mh,a,rat)
!c      CALL FLCLS2(NFL)
!c  999 CONTINUE
return
!xx    1 FORMAT(16I5)
!xx    2 FORMAT(4D20.10)
!xx  130 FORMAT(1H ,'PROGRAM 74.6.2. BISPEC / BISPECTRUM COMPUTATION.')
!xx  180 FORMAT(//1H ,4X,'C(I,J): THIRD ORDER MOMENTS'/)
!xx  129 FORMAT(1H ,'N=',I5,5X,'MH=',I5)
!xx  160 FORMAT(/1H ,'REAL PART OF BISPECTRUM B(I,J)'/)
!xx  161 FORMAT(/1H ,'IMAGINARY PART OF BISPECTRUM B(I,J)'/)
!xx  163 FORMAT(1H ,'CC(I):',1X,'AUTOCOVARIANCES')
!xx  165 FORMAT(1H ,I5,D16.5)
end
!
subroutine subca(cl,ca,mh,isw)
  use timsac_kinds, only: dp
  implicit none
!     ROWWISE SMOOTHING
!     CA(I,J)=
!              (CL(I-1,J)+2.0*CL(I,J)+CL(I+1,J))/4.0
!     ISW=0: COS
!     ISW=1: SIN
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION CL(51,51),CA(51,51)
!xx      DIMENSION CL(MH+1,MH+1),CA(MH+1,MH+1)
integer mh, isw
real(dp) cl(mh+1,mh+1), ca(mh+1,mh+1)
! local
integer i, i3, ii, ij, il, j, j1, j2, jj, jjl, l1
real(dp) cst2, cst4, cst4i
cst2=2.0d-00
cst4=4.0d-00
cst4i=1.0d-00/cst4
l1=mh+1
!     ON AND ABOVE X-AXIS COMPUTATION
jjl=mh/2+1
do 10 jj=1,jjl
j=jj-1
if(j.le.1) go to 11
ij=j
go to 12
11 ij=2
12 il=mh-j
do 20 ii=ij,il
ca(ii,jj)=&
&(cl(ii-1,jj)+cst2*cl(ii,jj)+cl(ii+1,jj))*cst4i
i=ii-1
20 continue
10 continue
!     BELOW X-AXIS
do 40 j1=1,2
jj=mh/2+j1+1
j2=j1+2
do 41 ii=j2,mh
i3=ii-j1
if(isw.eq.1) go to 42
ca(ii,jj)=ca(i3,j1+1)
i=ii-1
j=-j1
go to 41
42 ca(ii,jj)=-ca(i3,j1+1)
i=ii-1
j=-j1
41 continue
40 continue
return
end
!
subroutine subcb(ca,cb,mh)
  use timsac_kinds, only: dp
  implicit none
!     COLUMNWISE SMOOTHING
!     CB(I,J)=
!             (CA(I,J-1)+2.0*CA(I,J)+CA(I,J+1))/4.0
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION CA(51,51),CB(51,51)
!xx      DIMENSION CA(MH+1,MH+1),CB(MH+1,MH+1)
integer mh
real(dp) ca(mh+1,mh+1), cb(mh+1,mh+1)
! local
integer i, im1, j, jc, jj, jj1, jm1, l1, mh1, mhi, mhj
real(dp) cst2, cst4, cst4i
!     ON AND ABOVE 1-AXIS
cst2=2.0d-00
cst4=4.0d-00
cst4i=1.0d-00/cst4
l1=mh+1
mhj=mh/2
do 10 j=2,mhj
mhi=mh-j
do 11 i=j,mhi
cb(i,j)=&
&(ca(i,j-1)+cst2*ca(i,j)+ca(i,j+1))*cst4i
im1=i-1
jm1=j-1
11 continue
10 continue
!     ON 0-AXIS
mh1=mh-1
jj=mh/2+2
do 12 i=3,mh1
cb(i,1)=&
&(ca(i,jj)+cst2*ca(i,1)+ca(i,2))*cst4i
im1=i-1
jc=0
12 continue
!     ON (-1)-AXIS
jj1=mhj+1
do 14 i=4,mh
cb(i,jj1)=&
&(ca(i,jj+1)+cst2*ca(i,jj)+ca(i,1))*cst4i
im1=i-1
jc=-1
14 continue
return
end
!
!c      SUBROUTINE SUBCD(CB,CD,MH)
subroutine subcd(cb,cd,mh,b)
  use timsac_kinds, only: dp
  implicit none
!     OBLIQUE SMOOTHING
!     CD(I,J)=
!             (CB(I-1,J-1)+2.0*CB(I,J)+CB(I+1,J+1))/4.0
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      INTEGER KANA1 / ' B  ' /
!c      DIMENSION CB(51,51),CD(51,51)
!xx      DIMENSION CB(MH+1,MH+1),CD(MH+1,MH+1)
!xx      DIMENSION B(0:MH,0:MH)
integer mh
real(dp) cb(mh+1,mh+1), cd(mh+1,mh+1), b(0:mh,0:mh)
! local
integer i, il, im1, j, jc, jj1, jl, jm1, l1, mh3, mh4
real(dp) cst2, cst4, cst4i, cz
cst2=2.0d-00
cst4=4.0d-00
cst4i=1.0d-00/cst4
cz=0.0d0
!c      CALL BISPPR (KANA1,0,CZ,0)
l1=mh+1
!     ON AND ABOVE 2-AXIS
jl=mh/2-1
do 10 j=3,jl
il=mh-2-j
do 11 i=j,il
cd(i,j)=&
&(cb(i-1,j-1)+cst2*cb(i,j)+cb(i+1,j+1))*cst4i
im1=i-1
jm1=j-1
cz=cd(i,j)
!c      CALL BISPPR (IM1,JM1,CZ,1)
b(im1,jm1)=cz
11 continue
10 continue
!     ON 1-AXIS
mh4=mh-4
do 12 i=4,mh4
cd(i,2)=&
&(cb(i-1,1)+cst2*cb(i,2)+cb(i+1,3))*cst4i
im1=i-1
jc=1
cz=cd(i,2)
!c      CALL BISPPR(IM1,JC,CZ,1)
b(im1,jc)=cz
12 continue
!     ON 0-AXIS
mh3=mh-3
jj1=mh/2+1
do 14 i=5,mh3
cd(i,1)=&
&(cb(i-1,jj1)+cst2*cb(i,1)+cb(i+1,2))*cst4i
im1=i-1
jc=0
cz=cd(i,1)
!c      CALL BISPPR (IM1,JC,CZ,1)
b(im1,jc)=cz
14 continue
return
end
!
!c      SUBROUTINE SUBQ1(CL,SL,P,N,MH)
subroutine subq1(cl,sl,p,n,mh,a,rat)
  use timsac_kinds, only: dp
  implicit none
!     A MEASURE OF COHERENCE BETWEEN X(I)X(J) AND X(I+J) IS GIVEN BY Q1(
!     WHICH IS DEFINED BY
!     Q1(I,J)=(CL(I,J)**2+SL(I,J)**2)/(P(I)P(J)P(I+J)MH).
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      INTEGER KANA2 / 'Q1  ' /
!c      DIMENSION CL(51,51),SL(51,51),P(51)
!xx      DIMENSION CL(MH+1,MH+1),SL(MH+1,MH+1),P(MH+1)
!xx      DIMENSION A(0:MH,0:MH)
integer n, mh
real(dp) cl(mh+1,mh+1), sl(mh+1,mh+1), p(mh+1),&
&a(0:mh,0:mh), rat
! local
integer i, i4, il, im1, j, jc, jj1, jl, jm1, l1, mh3, mh4
real(dp) h, cz, an, as3, cst075
!
l1=mh+1
h = mh
cz =0.0d-00
!c      WRITE(6,185)
!c      CALL BISPPR (KANA2,0,CZ,0)
!     ON AND ABOVE 2-AXIS
jl=mh/2-1
!xx      DO 10 J=3,JL
do 11 j=3,jl
il=mh-2-j
do 10 i=j,il
i4=i+j-1
cl(i,j)=(cl(i,j)*cl(i,j)+sl(i,j)*sl(i,j))/p(i)/p(j)/p(i4)/h
im1=i-1
jm1=j-1
cz=cl(i,j)
!c      CALL BISPPR (IM1,JM1,CZ,1)
a(im1,jm1)=cz
10 continue
11 continue
!     ON 1-AXIS
mh4=mh-4
do 12 i=4,mh4
i4=i+1
cl(i,2)=(cl(i,2)*cl(i,2)+sl(i,2)*sl(i,2))/p(i)/p(2)/p(i4)/h
im1=i-1
jc=1
cz=cl(i,2)
!c      CALL BISPPR (IM1,JC,CZ,1)
a(im1,jc)=cz
12 continue
!     ON 0-AXIS
mh3=mh-3
jj1=mh/2+1
do 14 i=5,mh3
im1=i-1
cl(i,1)=(cl(i,1)*cl(i,1)+sl(i,1)*sl(i,1))/p(i)/p(i)/p(1)/h
jc=0
cz=cl(i,1)
!c      CALL BISPPR (IM1,JC,CZ,1)
a(im1,jc)=cz
14 continue
an=n
as3 =3.0d-00
cst075=0.75d-00
rat=(h/an)*cst075*cst075/dsqrt(as3)
!c      WRITE(6,186) RAT
return
!xx  181 FORMAT(1H ,2I5,D17.5)
!xx  185 FORMAT(/1H ,'Q1(I,J)=(CL(I,J)**2+SL(I,J)**2)/(P(I)P(J)P(I+J)MH);',
!xx     1'A MESURE OF COHERENCE BETWEEN X(I)*X(J) AND X(I+J)'/)
!xx  186 FORMAT(/1H ,'APPROXIMATE EXPECTED VALUE OF Q1(I,J) UNDER ',
!xx     A'GAUSSIAN ASSUMPTION;',D17.5,'.')
end
!
!c      SUBROUTINE SAUSP1(CXX,P1,N,LAGH3,LAGH1)
subroutine sausp1(cxx,p1,p2,q,n,lagh3,lagh1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES POWER SPECTRUM.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION CXX(51),FC(51),P1(51),P2(51),Q(51),P3(51)
!c      DIMENSION A1(10),A2(10)
!xx      DIMENSION CXX(LAGH1),FC(LAGH1),P1(LAGH1),P2(LAGH1),Q(LAGH1)
!xx      DIMENSION A1(2),A2(3)
integer n, lagh3, lagh1
real(dp) cxx(lagh1), p1(lagh1), p2(lagh1), q(lagh1)
! local
integer i, lagh, lagh0, mla1, mla2
real(dp) fc(lagh1), a1(2), a2(3)
!     WINDOW W1 DEFINITION
mla1=2
a1(1)=0.5d-00
a1(2) =0.25d-00
!     WINDOW W2 DEFINITION
mla2=3
a2(1)=0.625d-00
a2(2)=0.25d-00
a2(3)=-0.0625d-00
lagh=lagh1-1
lagh0=lagh3-1
do 10 i=2,lagh
!xx   10 CXX(I)=CXX(I)+CXX(I)
cxx(i)=cxx(i)+cxx(i)
10 continue
!     F-COS TRANSFORMATION
!     COMMON SUBROUTINE CALL
call fgerco(cxx,lagh1,fc,lagh1)
!     SPECTRUM SMOOTHING BY WINDOW W1
!     COMMON SUBROUTINE CALL
call ausp(fc,p1,lagh1,a1,mla1)
!     SPECTRUM SMOOTHING BY WINDOW W2
!     COMMON SUBROUTINE CALL
call ausp(fc,p2,lagh1,a2,mla2)
!     TEST STATISTICS COMPUTATION
!     COMMON SUBROUTINE CALL
call signif(p1,p2,q,lagh1,n)
!     AUTO SPECTRUM AND TEST STATISTICS PRINT OUT
!c      WRITE(6,66) N,LAGH
!c      WRITE(6,63) A1(1),A1(2)
!c      WRITE(6,64) A2(1),A2(2),A2(3)
!c      WRITE(6,67)
!     COMMON SUBROUTINE CALL
!c      CALL PRCOL3(P1,P2,Q,1,LAGH1,1)
return
!xx   63 FORMAT(1H ,'WINDOW W1',5X,'A(0)=',F10.5,5X,'A(1)=',F10.5)
!xx   64 FORMAT(1H ,'WINDOW W2',5X,'A(0)=',F10.5,5X,'A(1)=',F10.5,5X,'A(2)=
!xx     1',F10.5)
!xx   66 FORMAT(//1H ,14HPOWER SPECTRUM,5X,2HN=,I5,5X,5HLAGH=,I5)
!xx   67 FORMAT(/1H ,4X,1HI,8X,8HPOWER W1,6X,8HPOWER W2,2X,12HSIGNIFICANCE)
end
