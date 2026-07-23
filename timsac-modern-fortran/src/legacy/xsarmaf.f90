! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine xsarmaf( ys,n,iq,ip,p01,g1,tl1,p02,g2,alphb,alpha,tl2,&
&sigma2 )
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM XSARMA                                                    
!.......................................................................
!.....PLANNED BY H.AKAIKE...............................................
!.....DESIGNED BY H.AKAIKE..............................................
!.....PROGRAMMED BY E.ARAHATA...........................................
!.....ADDRESS: THE INSTITUTE OF STATISTICAL MATHEMATICS, 4-6-7 MINAMI-AZ
!..............MINATO-KU, TOKYO 106, JAPAN..............................
!.....DATE OF THE LATEST REVISION:  MAR. 6,1979.........................
!.......................................................................
!.....THIS PROGRAM WAS ORIGINALLY PUBLISHED IN "TIMSAC-78", BY H.AKAIKE,
!.....G.KITAGAWA, E.ARAHATA AND F.TADA, COMPUTER SCIENCE MONOGRAPHS, NO.
!.....THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO, 1979.............
!.......................................................................
!     TIMSAC 78.5.2                                                     
!     __                                 _      __ __                   
!     EXACT MAXIMUM LIKELIHOOD METHOD OF SCALAR AR-MA MODEL FITTING     
!                                                                       
!-----------------------------------------------------------------------
!     THIS PROGRAM PRODUCES EXACT MAXIMUM LIKELIHOOD ESTIMATES OF THE   
!     PARAMETERS OF A SCALAR AR-MA MODEL.                               
!                                                                       
!     THE AR-MA MODEL IS GIVEN BY                                       
!                                                                       
!     Y(I)+B(1)Y(I-1)+...+B(IQ)Y(I-IQ)=X(I)+A(1)X(I-1)+...+A(IP)X(I-IP),
!                                                                       
!     WHERE X(I) IS A ZERO MEAN WHITE NOISE.                            
!                                                                       
!     ----------------------------------------------------------------- 
!       REFERENCE:                                                      
!          H.AKAIKE(1978), "COVARIANCE MATRIX COMPUTATION OF THE STATE  
!          VARIABLE OF A STATIONARY GAUSSIAN PROCESS,",  RESEARCH MEMO. 
!          NO.139, THE INSTITUTE OF STATISTICAL MATHEMATICS; TOKYO.     
!          TO BE PUBLISHED IN ANN. INST. STATIST. MATH..                
!     ----------------------------------------------------------------- 
!     THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:    
!             SDATPR                                                    
!             SMINOP                                                    
!             SUBRST                                                    
!     ----------------------------------------------------------------- 
!     THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE SDATPR:          
!          IQ:  AR-ORDER                                                
!          (B(I),I=1,IQ):  INITIAL ESTIMATES OF AR COEFFICIENTS         
!          IP:  MA-ORDER                                                
!          (A(I),I=1,IP):  INITIAL ESTIMATES OF MA COEFFICIENTS         
!          TITL:  TITLE OF THE DATA                                     
!          N:  DATA LENGTH                                              
!          (DFORM(I),I=1,20):  INPUT FORMAT SPECIFICATION IN ONE CARD,  
!                              EXAMPLE,                                 
!                              (8F10.4)                                 
!          (Y(I),I=1,N):  ORIGINAL DATA                                 
!                                                                       
!     OUTPUTS:                                                          
!          TITL:  TITLE OF THE DATA                                     
!          IQ:  AR-ORDER                                                
!          (B(I),I=1,IQ):  AR COEFFICIENTS                              
!          IP:  MA-ORDER                                                
!          (A(I),I=1,IP):  MA COEFFICIENTS                              
!          SIGMA2:  WHITE NOISE VARIANCE                                
!-----------------------------------------------------------------------
!
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      REAL*4 TITL(20)                                                   
!c      DIMENSION Y(2000),P0(50)                                          
!xx      DIMENSION YS(N), Y(N), P01(IP+IQ)
!xx      DIMENSION P02(IP+IQ), G1(IP+IQ), G2(IP+IQ)
!xx      DIMENSION ALPHB(IQ), ALPHA(IP)
integer n, iq, ip
real(dp) ys(n), p01(ip+iq), g1(ip+iq), tl1, p02(ip+iq),&
&g2(ip+iq), alphb(iq), alpha(ip), tl2, sigma2
! local
real(dp) y(n)
!
!c      CHARACTER(100)  IFLNAM,OFLNAM,MFLNAM
!c      CALL FLNAM3( IFLNAM,OFLNAM,MFLNAM,NFL )
!c      IF ( NFL.EQ.0 ) GO TO 999
!c      IF ( (NFL.EQ.2).OR.(NFL.EQ.4) ) THEN
!c         OPEN (6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR)
!c      ELSE
!c         CALL SETWND
!c      END IF
!c      OPEN( 5,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
!c      IF ((NFL.EQ.3) .OR. (NFL.EQ.4)) THEN
!c         OPEN (7,FILE=MFLNAM,ERR=920,IOSTAT=IVAR)
!c      ELSE
!c         OPEN (7,FILE='xsarma.out',ERR=930,IOSTAT=IVAR)
!c      END IF
!
!                                                                       
!-----  READ IN AND PRINT OUT OF INITIAL CONDITION  -----               
!c      CALL SDATPR(TITL,Y,N,P0,IQ,IP)                                    
!xx      CALL SDATPR(YS,Y,N,P01,IQ,IP)
call sdatpr(ys,y,n)
!c      CLOSE( 5 )
!c#ifdef __linux__
!cC     reopen #5 as stdin
!c      OPEN(5, FILE='/dev/fd/0')
!c#endif
!cC /* __linux */
!                                                                       
!-----  MINIMIZATION OF TL(=-2)LOG LIKELIHOOD)                          
!       BY DAVIDON'S VARIANCE ALGORITHM  -----                          
!c      CALL SMINOP(TL,SIGMA2,Y,N,P0,IQ,IP)                               
!x      IPRNT = 0
!     IPRNT=0:  NOT TO PRINT OUT INTERMEDIATE RESULTS
!     IPRNT=1:  TO PRINT OUT INTERMEDIATE RESULTS ( LU: UNIT NUMBER)
!x      CALL SMINOP( TL1,TL2,SIGMA2,Y,N,P01,G1,P02,G2,ALPHB,ALPHA,IQ,IP,
!x     *             IPRNT,LU )
call sminop( tl1,tl2,sigma2,y,n,p01,g1,p02,g2,alphb,alpha,iq,ip)
!                                                                       
!-----  PRINT AND PUNCH OUT OF FINAL RESULT  -----                      
!c      CALL SUBRST(TITL,TL,SIGMA2,P0,IQ,IP)                              
!c      GO TO 990
!                                                                       
!
!c  900 CONTINUE
!c      WRITE(6,690) IVAR,OFLNAM
!c  690 FORMAT(' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
!c      GO TO 999
!
!c  910 CONTINUE
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
!c#ifdef __linux__
!cC     reopen #6 as stdout
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) OPEN(6, FILE='/dev/fd/1')
!c#endif
!cC /* __linux__ */
!c      WRITE(6,691) IVAR,IFLNAM
!c  691 FORMAT(' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
!c      GO TO 999
!
!c  920 CONTINUE
!c      CLOSE(5)
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
!c#ifdef __linux__
!cC     reopen #5 as stdin, #6 as stdout
!c      OPEN(5, FILE='/dev/fd/0')
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) OPEN(6, FILE='/dev/fd/1')
!c#endif
!cC /* __linux__ */
!c      WRITE(6,692) IVAR,MFLNAM
!c  692 FORMAT(' !!! Intermediate_Data_File OPEN ERROR ',I8,//,5X,100A)
!c      GO TO 999
!
!c  930 CONTINUE
!c      CLOSE(5)
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
!c#ifdef __linux__
!cC     reopen #5 as stdin, #6 as stdout
!c      OPEN(5, FILE='/dev/fd/0')
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) OPEN(6, FILE='/dev/fd/1')
!c#endif
!cC /* __linux__ */
!c      WRITE(6,693) IVAR
!c  693 FORMAT(' !!! xsarma.out  OPEN ERROR ',I8)
!c      GO TO 999
!
!c  990 CONTINUE
!c      CLOSE( 7 )
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
!c#ifdef __linux__
!cC     reopen #6 as stdout
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) OPEN(6, FILE='/dev/fd/1')
!c#endif
!cC /* __linux__ */
!c  999 CONTINUE
return
end
!c      SUBROUTINE ARCHCK(A,M,ICOND)                                      
subroutine archck(a,alph,m,icond)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!-----------------------------------------------------------------------
!     THIS SUBROUTINE CHECKS STABILITY OF AR OR MA PART.                
!                                                                       
!     INPUTS:                                                           
!          (A(I),I=1,M):  AR OR MA COEFFICIENTS                         
!          M:  AR- OR MA-ORDER                                          
!                                                                       
!     OUTPUTS:                                                          
!          (A(I),I=1,M):  AR OR MA COEFFICIENTS                         
!          ICOND:  ICOND=0, WHEN STABLE                                 
!                  ICOND=1, WHEN NOT STABLE                             
!-----------------------------------------------------------------------
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      DIMENSION A(50),B(50),ALPH(50)                                    
!c      COMMON /CMALPH/ALPH                                               
!xx      DIMENSION A(M),B(M),ALPH(M)                                    
integer m, icond
real(dp) a(m), alph(m)
! local
integer i, ik, im1, ip1, j, k, ki, mi, mip1
real(dp) b(m), cst0, cst1, cst099, at, ct, d
data cst0,cst1,cst099/0.0d-00,1.0d-00,0.99999d-00/
!                                                                       
!                                                                       
!                                                                       
!-----  STABILITY CHECK  -----                                          
!c      DO 10 I=1,50                                                      
!xx      DO 10 I=1,M
!xx   10 B(I)=CST0                                                         
b(1:m)=cst0
do 20 i=1,m
mi=m-i
mip1=mi+1
at=a(mip1)
if(dabs(at).lt.cst099) go to 210
icond=1
at=cst099*at/dabs(at)
210 alph(mip1)=at
if(mi.eq.0) go to 20
ct=cst1/(cst1-at**2)
do 25 j=1,mi
ki=mip1-j
b(j)=a(ki)
25 continue
do 30 k=1,mi
a(k)=(a(k)-at*b(k))*ct
30 continue
20 continue
!                                                                       
!-----  RECOVERY OF COEFFICIENTS  -----                                 
do 400 i=1,m
d=alph(i)
a(i)=d
if(i.eq.1) go to 410
im1=i-1
do 420 j=1,im1
a(j)=a(j)+d*b(j)
420 continue
410 continue
ip1=i+1
do 421 k=1,i
ik=ip1-k
!xx  421 B(K)=A(IK)
b(k)=a(ik)
421 continue
400 continue
!                                                                       
return
end
!c      SUBROUTINE  FUNCT2( F,SD,Y,N,P0,IQ,IP )                           
subroutine  funct2( f,sd,y,n,p0,iq,ip,ir )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!-----------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES SD AND F(=(-2)LOG LIKELIHOOD) BY         
!     A PROCEDURE OF MORF, SIDHU AND KAILATH (IEEE TRANS. AUTOMAT.      
!     CONTR., AC19, 315-323,1974).                                      
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
!             SUBPM                                                     
!       ----------------------------------------------------------------
!                                                                       
!     INPUTS:                                                           
!          N:  DATA LENGTH                                              
!          (Y(I),I=1,N):  ORIGINAL DATA, MEAN DELETED                   
!          IQ:  AR-ORDER                                                
!          IP:  MA-ORDER                                                
!          (P0(I),I=1,IPQ):  THE VECTOR OF AR AND MA COEFFICIENTS (IPQ=I
!                                                                       
!     OUTPUTS:                                                          
!          F:  (-2)LOG LIKELIHOOD                                       
!          SD:  WHITE NOISE VARIANCE                                    
!-----------------------------------------------------------------------
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      DIMENSION Y(2000)                                                 
!c      DIMENSION A(50),B(50),P(51,51),P0(50)                             
!c      DIMENSION P1(51),AKI(51),YY(51),PY(51),Z(51),PZ(51)               
!xx      DIMENSION Y(N)                                                 
!xx      DIMENSION A(IR),B(IR),P(IR,IR),P0(IP+IQ)                             
!xx      DIMENSION P1(IR),AKI(IR),YY(IR),PY(IR),Z(IR),PZ(IR)
integer n,  iq, ip, ir
real(dp) f ,sd, y(n), p0(ip+iq)
! local
integer i, i1, ii, indx, indx1, irm1, irp1, ns
real(dp) a(ir), b(ir), p(ir,ir), p1(ir), aki(ir), yy(ir),&
&py(ir), z(ir), pz(ir), cst0, cst1, csta, rem,&
&sum, ami, em, seri, slr, y2, rew, eri, amy,&
&amri, yri, ew, an
!c      DATA A,B,AKI/151*0.0D-00/                                         
data cst0,cst1,csta/0.0d-00,1.0d-00,0.1d-05/
!                                                                       
!xx      DO 600 I=1,IR
!xx         A(I)=CST0
!xx         B(I)=CST0
!xx         AKI(I)=CST0
!xx  600 CONTINUE
a(1:ir)=cst0
b(1:ir)=cst0
aki(1:ir)=cst0
!                                                                       
if(iq.eq.0) go to 620
do 610 i=1,iq
b(i)=p0(i)
610 continue
620 if(ip.eq.0) go to 640
do 630 i=1,ip
ii=iq+i
a(i)=p0(ii)
630 continue
640 continue
!                                                                       
!c      IP1=IP+1                                                          
!c      IR=MAX0(IQ,IP1)                                                   
!                                                                       
!-----  P MATRIX COMPUTATION  -----                                     
call subpm(p,b,a,iq,ip,ir)
!                                                                       
!-----  INITIAL CONDITION  -----                                        
rem=p(1,1)
do 20 i=1,ir
p1(i)=p(i,1)
20 continue
irm1=ir-1
irp1=ir+1
if(irm1.le.0) go to 510
do 50 i=1,irm1
aki(i)=p1(i+1)
50 continue
510 continue
sum=cst0
do 60 i=1,ir
i1=irp1-i
sum=sum+b(i)*p1(i1)
60 continue
aki(ir)=-sum
ami=-cst1/rem
do 100 i=1,ir
yy(i)=aki(i)
100 continue
em=y(1)
do 110 i=1,ir
z(i)=cst0
110 continue
!                                                                       
seri=em**2/rem
slr=dlog(rem)
!******************************                                         
indx = 1
do 1000 ns=2,n
y2=yy(1)*yy(1)
!-----  RE COMPUTATION  -----                                           
rew=rem+ami*y2
indx=ns
!                                                                       
!-----  PHI*Z COMPUTATION  -----                                        
if(irm1.le.0) go to 520
do 210 i=1,irm1
pz(i)=z(i+1)
210 continue
520 continue
sum=cst0
do 220 i=1,ir
i1=irp1-i
sum=sum+b(i)*z(i1)
220 continue
pz(ir)=-sum
!----- Z COMPUTATION  -----                                             
eri=em/rem
do 230 i=1,ir
z(i)=pz(i)+eri*aki(i)
230 continue
!                                                                       
!-----  PHI*YY COMPUTATION  -----                                       
if(irm1.le.0) go to 530
do 310 i=1,irm1
py(i)=yy(i+1)
310 continue
530 continue
sum=cst0
do 320 i=1,ir
i1=irp1-i
sum=sum+b(i)*yy(i1)
320 continue
py(ir)=-sum
!                                                                       
!-----  K COMPUTATION  -----                                            
amy=ami*yy(1)
do 330 i=1,ir
aki(i)=aki(i)+amy*py(i)
330 continue
!                                                                       
!-----  M COMPUTATION  -----                                            
amri=ami/rem
ami=ami*(cst1+amri*y2)
!                                                                       
!-----  YY COMPUTATION  -----                                           
yri=yy(1)/rew
do 400 i=1,ir
yy(i)=py(i)-yri*aki(i)
400 continue
!                                                                       
!-----  E COMPUTATION  -----                                            
ew=y(ns)-z(1)
!                                                                       
!-----  SERI, SLR COMPUTATION  -----                                    
seri=seri+ew**2/rew
slr=slr+dlog(rew)
rem=rew
em=ew
if(dabs(rew-cst1).lt.csta) go to 1100
1000 continue
!                                                                       
1100 continue
if(indx.ge.n) go to 1500
indx1=indx+1
do 1110 ns=indx1,n
!-----  PHI*Z  -----                                                    
if(irm1.le.0) go to 540
do 1210 i=1,irm1
pz(i)=z(i+1)
1210 continue
540 continue
sum=cst0
do 1220 i=1,ir
i1=irp1-i
sum=sum+b(i)*z(i1)
1220 continue
pz(ir)=-sum
!                                                                       
!-----  Z COMPUTATION  -----                                            
do 1230 i=1,ir
z(i)=pz(i)+em*aki(i)
1230 continue
!                                                                       
!-----  E COMPUTATION  -----                                            
ew=y(ns)-z(1)
!                                                                       
!-----  SERI COMPUTATION  -----                                         
seri=seri+ew**2
em=ew
1110 continue
1500 continue
!******************************                                         
!                                                                       
!                                                                       
an=n
sd=seri/an
f=slr+an*dlog(sd)
!                                                                       
if(iq.eq.0) go to 1620
do 1610 i=1,iq
p0(i)=b(i)
1610 continue
1620 if(ip.eq.0) go to 1640
do 1630 i=1,ip
ii=iq+i
p0(ii)=a(i)
1630 continue
1640 continue
!xx 2100 CONTINUE                                                          
!                                                                       
return
end
!c      SUBROUTINE MSDAV2(PHAI,SIGMA2,G,C,Y,N,X,IQ,IP,ISWRO,IPRNT)        
!x      SUBROUTINE MSDAV2(PHAI,SIGMA2,G,C,Y,N,X,IQ,IP,ISWRO,VD,IPRNT,LU)
subroutine msdav2(phai,sigma2,g,c,y,n,x,iq,ip,iswro,vd)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!-----------------------------------------------------------------------
!     DAVIDON'S (MINIMIZATION) PROCEDURE                                
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             SGRAD                                                     
!             ARCHCK                                                    
!                                                                       
!     INPUTS:                                                           
!          N:  DATA LENGTH                                              
!          (Y(I),I=1,N):  ORIGINAL DATA, MEAN DELETED                   
!          IQ:  AR-ORDER                                                
!          IP:  MA-ORDER                                                
!          (X(I),I=1,IPQ):  THE VECTOR OF AR AND MA COEFFICIENTS (IPQ=IP
!          PHAI:  (-2)LOG LIKELIHOOD                                    
!          SIGMA2:  WHITE NOISE VARIANCE                                
!          (G(I),I=1,IPQ):  GRADIENT                                    
!          (C(I),I=1,IPQ):  CORRECTION TERM                             
!          ISWRO:  ITERATION COUNT OF SUBROUTINE MSDAV2                 
!          ((VD(I,J),I=1,IPQ),J=1,IPQ):  INVERSE OF HESSIAN             
!                                                                       
!     OUTPUTS:                                                          
!          PHAI:  NEW PHAI                                              
!          SIGMA2:  NEW SIGMA2                                          
!          (G(I),I=1,IPQ):  NEW GRADIENT                                
!          (X(I),I=1,IPQ): NEW X(I)                                     
!          ISWRO:  NEW ISWRO                                            
!          ((VD(I,J),I=1,IPQ),J=1,IPQ):  NEW INVERSE OF HESSIAN         
!                                                                       
!     PARAMETERS:                                                       
!          IPRNT=0:  NOT TO PRINT OUT INTERMEDIATE RESULTS              
!          IPRNT=1:  TO PRINT OUT INTERMEDIATE RESULTS                  
!-----------------------------------------------------------------------
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!xx      REAL*8 MAXVD                                                      
!c      DIMENSION VD(50,50),X(50),G(50),SX(50),SG(50),SR(50),C(50)        
!c      DIMENSION Y(2000),SSX(50)                                         
!c      COMMON /COM50/VD                                                  
!xx      DIMENSION G(IP+IQ),C(IP+IQ),Y(N),X(IP+IQ),VD(IP+IQ,IP+IQ)
!xx      DIMENSION SX(IP+IQ),SG(IP+IQ),SR(IP+IQ),SSX(IP+IQ)                                         
!xx      DIMENSION ALPH(IP+IQ)
integer n, iq, ip, iswro
real(dp) phai, sigma2, g(ip+iq), c(ip+iq), y(n), x(ip+iq),&
&vd(ip+iq,ip+iq)
! local
integer i, icond, ii, iphai, ipq, ipq2, iram, isphai, itn, itns,&
&j
real(dp) maxvd, sx(ip+iq), sg(ip+iq), sr(ip+iq),&
&ssx(ip+iq), alph(ip+iq), cst0, cst1, cst4, cst10,&
&consta, constb, eps3, eps4, vdn, sci, sum,&
&sro, gsr, dgam, dgam1, ram, ramsro, ramt,&
&sphai, ram1, consdr, sd
data cst0,cst1,cst4,cst10/0.0d-00,1.0d-00,4.0d-00,10.0d-00/
data consta,constb,eps3,eps4/0.5d-00,2.0d-00,0.1d-05,0.1d-05/
!c      DATA SSX/50*0.0D-00/                                              
!                                                                       
!                                                                       
isphai=0
iphai=1
ipq=iq+ip
ipq2=ipq+ipq
!                                                                       
itns=0
150 continue
itns=itns+1
itn=0
!                                                                       
1210 continue
icond=0
!-----  SX=X-C  -----                                                   
do 210 i=1,ipq
!xx  210 SX(I)=X(I)-C(I)
sx(i)=x(i)-c(i)
210 continue
!                                                                       
!x      IF(IPRNT.EQ.0) GO TO 3200                                         
!c      WRITE(6,3900)                                                     
!c      WRITE(6,3910) (C(I),I=1,IPQ)                                      
!x      WRITE(LU,3900)                                                     
!x      WRITE(LU,3910) (C(I),I=1,IPQ)                                      
!x 3200 CONTINUE                                                          
!                                                                       
!                                                                       
!-----  GRADIENT COMPUTATION  -----                                     
!                                                                       
!xx 4000 CONTINUE                                                          
icond=0
if(iq.le.0) go to 4510
do 4500 i=1,iq
!xx 4500 SSX(I)=SX(I)
ssx(i)=sx(i)
4500 continue
!c      CALL ARCHCK(SSX,IQ,ICOND)                                         
call archck(ssx,alph,iq,icond)
4510 if(ip.le.0) go to 4600
do 4520 i=1,ip
ii=iq+i
!xx 4520 SSX(I)=SX(II)
ssx(i)=sx(ii)
4520 continue
!c      CALL ARCHCK(SSX,IP,ICOND)                                         
call archck(ssx,alph,ip,icond)
4600 continue
!                                                                       
if(icond.eq.0) go to 309
!                                                                       
!x      IF(IPRNT.EQ.0) GO TO 1220                                         
!c      WRITE(6,4700)                                                     
!x      WRITE(LU,4700)                                                     
1220 continue
!                                                                       
itn=itn+1
!                                                                       
!------------------------------                                         
maxvd=cst0
do 4900 i=1,ipq
if(vd(i,i).gt.maxvd)maxvd=vd(i,i)
4900 continue
vdn=maxvd/cst4
do 304 i=1,ipq
do 303 j=1,ipq
!xx  303 VD(I,J)=VD(I,J)/CST10
vd(i,j)=vd(i,j)/cst10
303 continue
!xx  304 VD(I,I)=VD(I,I)+VDN
vd(i,i)=vd(i,i)+vdn
304 continue
do 306 i=1,ipq
sci=cst0
do 305 j=1,ipq
!xx  305 SCI=SCI+VD(I,J)*G(J)
sci=sci+vd(i,j)*g(j)
305 continue
!xx  306 C(I)=SCI
c(i)=sci
306 continue
!------------------------------                                         
!                                                                       
go to 1210
309 continue
!c      CALL SGRAD(SPHAI,SD,SG,Y,N,SX,IQ,IP,IPRNT)                        
!x      CALL SGRAD(SPHAI,SD,SG,Y,N,SX,IQ,IP,IPRNT,LU)                        
call sgrad(sphai,sd,sg,y,n,sx,iq,ip)
if(icond.eq.1) go to 1220
if(itn.ge.10) go to 312
312 continue
!                                                                       
!-----  SR=V*SG  -----                                                  
do 310 i=1,ipq
sum=cst0
do 311 j=1,ipq
!xx  311 SUM=SUM+VD(I,J)*SG(J)                                             
sum=sum+vd(i,j)*sg(j)
311 continue
!xx  310 SR(I)=SUM
sr(i)=sum
310 continue
!                                                                       
!-----  SRO=(SG)'*(SR)  -----                                           
sro=0.0d-00
do 1050 i=1,ipq
!xx 1050 SRO=SRO+SG(I)*SR(I)
sro=sro+sg(i)*sr(i)
1050 continue
!                                                                       
!-----  DGAM=-G'*(SR)/SRO  -----                                        
gsr=0.0d-00
do 1060 i=1,ipq
!xx 1060 GSR=GSR+G(I)*SR(I)
gsr=gsr+g(i)*sr(i)
1060 continue
dgam=-gsr/sro
dgam1=dgam+cst1
dgam1=dabs(dgam1)+0.1d-70
ram=dabs(dgam)/dgam1
!                                                                       
if(ram.gt.consta) go to 430
ram=consta
iram=1
go to 470
!                                                                       
430 if(ram.lt.constb) go to 450
ram=constb
iram=-1
go to 470
!                                                                       
450 continue
iram=0
470 ramsro=(ram-cst1)/sro
!xx      DO 480 I=1,IPQ                                                    
do 481 i=1,ipq
ramt=ramsro*sr(i)
do 480 j=1,ipq
!xx  480 VD(I,J)=VD(I,J)+RAMT*SR(J)
vd(i,j)=vd(i,j)+ramt*sr(j)
480 continue
481 continue
if(phai.ge.sphai) go to 540
!                                                                       
!-----  SPHAI.GT.PHAI: TEST OF CORRECTION  -----                        
ram1=ram-cst1
if(dabs(ram1).lt.eps3) go to 555
consdr=dgam*ram1
do 550 i=1,ipq
!xx  550 C(I)=C(I)-CONSDR*SR(I)                                            
c(i)=c(i)-consdr*sr(i)
550 continue
iphai=0
if(sro.gt.eps4) go to 900
!     END OF ITERATION                                                  
555 iswro=iswro+1
go to 1000
!                                                                       
!-----  SHPAI.LE.PHAI: SUCCESSFUL REDUCTION  -----                      
540 do 560 i=1,ipq
x(i)=sx(i)
g(i)=sg(i)
!xx  560 C(I)=RAM*SR(I)                                                    
c(i)=ram*sr(i)
560 continue
phai=sphai
sigma2=sd
!                                                                       
!x      IF(IPRNT.EQ.0) GO TO 571                                          
!c      WRITE(6,570) PHAI                                                 
!x      WRITE(LU,570) PHAI                                                 
!x  571 CONTINUE                                                          
!                                                                       
iphai=1
!xx  800 CONTINUE
if(iram.ne.0) go to 901
if(sro.lt.eps4) go to 555
!     ITERATION CHECK                                                   
900 continue
isphai=(isphai+(1-iphai))*(1-iphai)
!                                                                       
if(isphai.gt.ipq2) go to 555
!                                                                       
go to 150
901 if(sro.lt.eps4) go to 555
go to 900
!     END OF MINIMIZATION                                               
1000 continue
!                                                                       
!xx 1001 RETURN                                                            
!xx  570 FORMAT(1H ,'NEW PHAI=',D20.10)                                    
!xx 3900 FORMAT(//1H ,'C(I)')                                              
!xx 3910 FORMAT(1H ,5D20.10)                                               
!xx 4700 FORMAT(1H ,'ON THE BOUNDARY')                                     
return
end
!c      SUBROUTINE SDATPR(TITL,Y,N,P0,IQ,IP)                              
!xx      SUBROUTINE SDATPR(YS,Y,N,P0,IQ,IP)                              
subroutine sdatpr(ys,y,n)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!-----------------------------------------------------------------------
!     THIS SUBROUTINE READS IN AND PRINTS OUT INITIAL CONDITION AND DELE
!     THE MEAN OF THE DATA.                                             
!                                                                       
!     THE FOLLOWING INPUTS ARE REQUIRED:                                
!          IQ:  AR-ORDER                                                
!          (B(I),I=1,IQ):  INITIAL ESTIMATES OF AR COEFFICIENTS         
!          IP:  MA-ORDER                                                
!          (A(I),I=1,IP):  INITIAL ESTIMATES OF MA COEFFICIENTS         
!          TITL:  TITLE OF THE DATA                                     
!          N:  DATA LENGTH                                              
!          (DFORM(I),I=1,20):  INPUT FORMAT SPECIFICATION IN ONE CARD,  
!                              EXAMPLE,                                 
!                              (8F10.4)                                 
!          (Y(I),I=1,N):  ORIGINAL DATA                                 
!                                                                       
!     THE AR-MA MODEL IS GIVEN BY                                       
!     Y(I)+B(1)Y(I-1)+...+B(IQ)Y(I-IQ)=X(I)+A(1)X(I-1)+...+A(IP)X(I-IP),
!     WHERE X(I) IS A ZERO MEAN WHITE NOISE.                            
!-----------------------------------------------------------------------
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      REAL*4 DFORM(20),TITL(20)                                         
!c      DIMENSION Y(2000),A(50),B(50),P0(50)                              
!c      DATA A,B/100*0.0D-00/                                             
!xx      DIMENSION YS(N), Y(N), P0(IP+IQ)
integer n
real(dp) ys(n), y(n)
! local
integer i
real(dp) cst0, an, sum, ymean
data cst0/0.0d-00/
!                                                                       
!-----  INITIAL CONDITION LOADING FOR AR-MA(IQ,IP)  -----               
!c      READ(5,1) IQ                                                      
!c      IF(IQ.LE.0) GO TO 4205                                            
!c      READ(5,2) (B(I),I=1,IQ)                                           
!c 4205 READ(5,1) IP                                                      
!c      IF(IP.LE.0) GO TO 4206                                            
!c      READ(5,2) (A(I),I=1,IP)                                           
!c 4206 CONTINUE                                                          
!                                                                       
!-----  P0 ARRANGEMENT  -----                                           
!c      IF(IQ.LE.0) GO TO 300                                             
!c      DO 200 I=1,IQ                                                     
!c      P0(I)=B(I)                                                        
!c  200 CONTINUE                                                          
!c  300 IF(IP.LE.0) GO TO 310                                             
!c      DO 210 I=1,IP                                                     
!c      II=IQ+I                                                           
!c      P0(II)=A(I)                                                       
!c  210 CONTINUE                                                          
!c  310 CONTINUE                                                          
!                                                                       
!-----  DATA INPUT  -----                                               
!c      READ(5,4) (TITL(I),I=1,20)                                        
!c      READ(5,1) N                                                       
!c      READ(5,4) (DFORM(I),I=1,20)                                       
!     ORIGINAL DATA INPUT AND PRINT OUT                                 
!c      READ(5,DFORM) (Y(I),I=1,N)                                        
do 320 i=1,n
!xx  320 Y(I)=YS(I)
y(i)=ys(i)
320 continue
!                                                                       
!c      WRITE(6,59)                                                       
!c      WRITE(6,162) (TITL(I),I=1,20)                                     
!c      WRITE(6,62) IQ                                                    
!c      IF(IQ.LE.0) GO TO 4215                                            
!c      WRITE(6,65) (B(I),I=1,IQ)                                         
!c 4215 CONTINUE                                                          
!c      WRITE(6,63) IP                                                    
!c      IF(IP.LE.0) GO TO 4216                                            
!c      WRITE(6,65) (A(I),I=1,IP)                                         
!c 4216 CONTINUE                                                          
!                                                                       
!c      WRITE(6,285)                                                      
!c      WRITE(6,164) N                                                    
!c      WRITE(6,610) (Y(I),I=1,N)                                         
!                                                                       
!-----  MEAN DELETION  -----                                            
an=n
sum=cst0
do 9 i=1,n
!xx    9 SUM=SUM+Y(I)                                                      
sum=sum+y(i)
9 continue
ymean=sum/an
do 10 i=1,n
!xx   10 Y(I)=Y(I)-YMEAN
y(i)=y(i)-ymean
10 continue
!                                                                       
return
!xx    1 FORMAT(16I5)                                                      
!xx    2 FORMAT(4D20.10)                                                   
!xx    4 FORMAT(20A4)                                                      
!xx   59 FORMAT( //,' PROGRAM TIMSAC 78.5.2',/,
!xx     *'   EXACT MAXIMUM LIKELIHOOD METHOD OF AR-MA MODEL FITTING;',
!xx     *'  SCALAR CASE',/,'   < AR-MA MODEL >',
!xx     *//,11X,'Y(I) + B(1)*Y(I-1) + ... + B(IQ)*Y(I-IQ)  =  ',
!xx     *'X(I) + A(1)*X(I-1) + ... + A(IP)*X(I-IP)',/,' WHERE',/,11X,
!xx     *'IQ:    AR-ORDER',/,11X,'IP:    MA-ORDER',/,11X,
!xx     *'X(I): ZERO MEAN WHITE NOISE' ) 
!xx   62 FORMAT(//1H ,5('-'),2X,'INITIAL AR(I)',2X,'IQ=',I3,2X,5('-'))     
!xx   63 FORMAT(//1H ,5('-'),2X,'INITIAL MA(I)',2X,'IP=',I3,2X,5('-'))     
!xx   65 FORMAT(1H ,5D20.10,/(1H ,5D20.10))                                
!xx  162 FORMAT( 1H ,' TITLE:',/,2X,20A4 )                                 
!xx  164 FORMAT(1H ,'N=',I5)                                               
!xx  285 FORMAT(///1H ,5('-'),2X,'ORIGINAL DATA',2X,5('-'))                
!xx  610 FORMAT(1H ,10F10.5,/(1H ,10F10.5))                                
end
!c      SUBROUTINE SGRAD(F0,SD,G,Y,N,P0,IQ,IP,IPRNT)                      
!x      SUBROUTINE SGRAD(F0,SD,G,Y,N,P0,IQ,IP,IPRNT,LU)                      
subroutine sgrad(f0,sd,g,y,n,p0,iq,ip)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!-----------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES AN APPROXIMATION TO GRADIENT BY DIFFERENC
!     THIS SUBROUTINE SHOULD EVENTUALLY BE REPLACED BY AN ANALYTIC EVALU
!     PROCEDURE OF GRADIENTS.                                           
!       ----------------------------------------------------------------
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             FUNCT2                                                    
!             ARCHCK                                                    
!       ----------------------------------------------------------------
!                                                                       
!     INPUTS:                                                           
!          N:  DATA LENGTH                                              
!          (Y(I),I=1,N):  ORIGINAL DATA, MEAN DELETED                   
!          IQ:  AR-ORDER                                                
!          IP:  MA-ORDER                                                
!          (P0(I),I=1,IPQ):  THE VECTOR OF AR AND MA COEFFICIENTS (IPQ=I
!                                                                       
!     OUTPUTS:                                                          
!          F0:  (-2)LOG LIKELIHOOD                                      
!          SD:  WHITE NOISE VARIANCE                                    
!          (G(I),I=1,IPQ):  GRADIENT                                    
!                                                                       
!     PARAMETERS:                                                       
!          EPSA:  ORDINATE DIFFERENCE FOR GRADIENT COMPUTATION BY DIFFER
!          IPRNT=0:  NOT TO PRINT OUT INTERMEDIATE RESULTS              
!          IPRNT=1:  TO PRINT OUT INTERMEDIATE RESULTS                  
!-----------------------------------------------------------------------
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      DIMENSION P0(50),P1(50),Y(2000),G(50),PP0(50)                     
!xx      DIMENSION P0(IP+IQ),P1(IP+IQ),Y(N),G(IP+IQ),PP0(IP+IQ)                     
!xx      DIMENSION ALPH(IP+IQ)
integer n, iq, ip
real(dp) f0, sd, g(ip+iq), y(n), p0(ip+iq)
! local
integer i, icond, ii, iii, ip1, ipq, ir, itr, j
real(dp) p1(ip+iq), pp0(ip+iq), alph(ip+iq), cst07,&
&epsa, epsas, f1, sdn
data epsa,cst07/0.1d-03,0.7d-00/
!                                                                       
!c      CALL FUNCT2(F0,SD,Y,N,P0,IQ,IP)                                   
ip1=ip+1
ir=max0(iq,ip1)
call funct2(f0,sd,y,n,p0,iq,ip,ir)
ipq=ip+iq
do 9 j=1,ipq
p1(j)=p0(j)
9 continue
!                                                                       
!-----  GRADIENT COMPUTATION  -----                                     
do 10 i=1,ipq
epsas=epsa
itr=1
4000 continue
icond=0
p1(i)=p0(i)+epsas
if(i.gt.iq) go to 4500
do 4010 ii=1,iq
pp0(ii)=p1(ii)
4010 continue
!c      CALL ARCHCK(PP0,IQ,ICOND)                                         
call archck(pp0,alph,iq,icond)
go to 5000
4500 continue
do 5010 ii=1,ip
iii=iq+ii
!xx 5010 PP0(II)=P1(III)                                                   
pp0(ii)=p1(iii)
5010 continue
!c      CALL ARCHCK(PP0,IP,ICOND)                                         
call archck(pp0,alph,ip,icond)
5000 continue
if(icond.eq.0) go to 4100
if(itr.lt.10) go to 4110
!     WRITE(6,4200)                                                     
return
4110 epsas=epsas*cst07
epsas=-epsas
itr=itr+1
go to 4000
4100 continue
!c      CALL FUNCT2(F1,SDN,Y,N,P1,IQ,IP)                                  
call funct2(f1,sdn,y,n,p1,iq,ip,ir)
g(i)=(f1-f0)/epsas
p1(i)=p0(i)
10 continue
!--------------------                                                   
!                                                                       
!x      IF(IPRNT.EQ.0) RETURN                                             
!c      WRITE(6,450)                                                      
!c      WRITE(6,520) (P0(I),I=1,IPQ)                                      
!c      WRITE(6,500)                                                      
!c      WRITE(6,520) (G(I),I=1,IPQ)                                       
!x      WRITE(LU,450)                                                      
!x      WRITE(LU,520) (P0(I),I=1,IPQ)                                      
!x      WRITE(LU,500)                                                      
!x      WRITE(LU,520) (G(I),I=1,IPQ)                                       
!                                                                       
return
!xx  450 FORMAT(1H ,'P0(I)')                                               
!xx  500 FORMAT(1H ,'GRADIENT G(I)')                                       
!xx  520 FORMAT(1H ,5D20.10)                                               
!xx 4200 FORMAT(1H ,'ICOND=1; GRADIENT UNOBTAINED')                        
end
!c      SUBROUTINE SMINOP(TL,SIGMA2,Y,N,P0,IQ,IP)                         
subroutine sminop( tl,tl2,sigma2,y,n,p0,g,p02,g2,alphb,alpha,iq,&
!x     *                   IP,IPRNT,LU )
&ip )
  use timsac_kinds, only: dp
  implicit none
! C                                                                       
!-----------------------------------------------------------------------
!     THIS SUBROUTINE CONTROLS THE MAXIMUM LIKELIHOOD COMPUTATION.      
!                                                                       
!       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
!             SGRAD                                                     
!             ARCHCK                                                    
!                                                                       
!     INPUTS:                                                           
!          N:  DATA LENGTH                                              
!          (Y(I),I=1,N):  ORIGINAL DATA, MEAN DELETED                   
!          IQ:  AR-ORDER                                                
!          IP:  MA-ORDER                                                
!          (P0(I),I=1,IPQ):  THE VECTOR OF AR AND MA COEFFICIENTS (IPQ=I
!                                                                       
!     OUTPUTS:                                                          
!          TL:  (-2)LOG LIKELIHOOD                                      
!          SIGMA2:  WHITE NOISE VARIANCE                                
!                                                                       
!     PARAMETERS:                                                       
!          IPRNT=0:  NOT TO PRINT OUT INTERMEDIATE RESULTS              
!          IPRNT=1:  TO PRINT OUT INTERMEDIATE RESULTS                  
!-----------------------------------------------------------------------
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!xx      REAL*8 MAXAB                                                      
!c      DIMENSION Y(2000),P0(50),G(50),HS(50,50),CR(50),PP0(50),ALPH(50)  
!xx      DIMENSION Y(N),P0(IP+IQ),P02(IP+IQ)
!xx      DIMENSION G(IP+IQ),G2(IP+IQ),HS(IP+IQ,IP+IQ),CR(IP+IQ)
!xx      DIMENSION PP0(IP+IQ),ALPH(IP+IQ),ALPHB(IQ),ALPHA(IP)
integer n, iq, ip
real(dp) tl, tl2, sigma2, y(n), p0(ip+iq), g(ip+iq),&
&p02(ip+iq), g2(ip+iq), alphb(iq), alpha(ip)
! local
integer i, icond, ii, ipq, iswro, j
real(dp) maxab, hs(ip+iq,ip+iq), cr(ip+iq), pp0(ip+iq),&
&alph(ip+iq), cst0, cst10, cst05, cst005, pab,&
&bn, sum
!c      COMMON /COM50/HS /CMALPH/ALPH                                     
!c      DATA PP0/50*0.0D-00/                                              
data cst0,cst10,cst05,cst005/0.0d-00,10.0d-00,0.1d-03,0.00005d-00/
!                                                                       
!                                                                       
!c      IPRNT=0                                                           
!                                                                       
ipq=ip+iq
!xx      DO 310 I=1,IPQ                                                    
!xx      G(I)=CST0                                                         
!xx      PP0(I)=CST0
!xx      DO 310 J=1,IPQ                                                    
!xx  310 HS(I,J)=CST0                                                      
g(1:ipq)=cst0
pp0(1:ipq)=cst0
hs(1:ipq,1:ipq)=cst0
!                                                                       
!                                                                       
!-----  INITIAL GRADIENT COMPUTATION  -----
!xx 4000 CONTINUE                                                          
icond=0
if(iq.le.0) go to 4510
do 4500 i=1,iq
!xx 4500 PP0(I)=P0(I)                                                      
pp0(i)=p0(i)
4500 continue
!c      CALL ARCHCK(PP0,IQ,ICOND)                                         
call archck(pp0,alph,iq,icond)
do 5000 i=1,iq
!xx 5000 P0(I)=PP0(I)                                                      
p0(i)=pp0(i)
5000 continue
4510 if(ip.le.0) go to 4800
do 4700 i=1,ip
ii=iq+i
!xx 4700 PP0(I)=P0(II)                                                     
pp0(i)=p0(ii)
4700 continue
!c      CALL ARCHCK(PP0,IP,ICOND)                                         
call archck(pp0,alph,ip,icond)
do 5100 i=1,ip
ii=iq+i
p0(ii)=pp0(i)
5100 continue
4800 continue
iswro=0
!c      CALL SGRAD(TL,SIGMA2,G,Y,N,P0,IQ,IP,IPRNT)                        
!x      CALL SGRAD(TL,SIGMA2,G,Y,N,P0,IQ,IP,IPRNT,LU)                        
call sgrad(tl,sigma2,g,y,n,p0,iq,ip)
!                                                                       
!c      WRITE(6,450)                                                      
!c      WRITE(6,520) (P0(I),I=1,IPQ)                                      
!c      WRITE(6,500)                                                      
!c      WRITE(6,520) (G(I),I=1,IPQ)                                       
!c      WRITE(6,550) TL                                                   
do 4850 i=1,ipq
p02(i)=p0(i)
g2(i)=g(i)
4850 continue
tl2=tl
!                                                                       
4890 continue
maxab=cst0
do 4900 i=1,ipq
!c      PAB=DABS(G(I))                                                    
pab=dabs(g2(i))
if(pab.gt.maxab) maxab=pab
4900 continue
!                                                                       
!                                                                       
!-----  INVERSE OF HESSIAN COMPUTATION  -----                           
bn=cst05/maxab
do 3010 i=1,ipq
do 3009 j=1,ipq
!xx 3009 HS(I,J)=HS(I,J)/CST10                                             
hs(i,j)=hs(i,j)/cst10
3009 continue
hs(i,i)=bn+hs(i,i)
3010 continue
!                                                                       
!x      IF(IPRNT.EQ.0) GO TO 3120                                         
!c      WRITE(6,3000)                                                     
!x      WRITE(LU,3000)                                                     
!x      DO 3100 I=1,IPQ                                                   
!c 3100 WRITE(6,3110) I,(HS(I,J),J=1,IPQ)                                 
!x 3100 WRITE(LU,3110) I,(HS(I,J),J=1,IPQ)                                 
!x 3120 CONTINUE                                                          
!                                                                       
!-----  CORRECTION TERM CR(X)=HS*G(X) COMPUTATION  -----                
do 900 i=1,ipq
sum=cst0
do 910 j=1,ipq
!c  910 SUM=SUM+HS(I,J)*G(J)                                              
!xx  910 SUM=SUM+HS(I,J)*G2(J)
sum=sum+hs(i,j)*g2(j)
910 continue
!xx  900 CR(I)=SUM                                                         
cr(i)=sum
900 continue
!                                                                       
!x      IF(IPRNT.EQ.0) GO TO 3920                                         
!c      WRITE(6,3900)                                                     
!c      WRITE(6,3910) (CR(I),I=1,IPQ)                                     
!x      WRITE(LU,3900)                                                     
!x      WRITE(LU,3910) (CR(I),I=1,IPQ)                                     
!x 3920 CONTINUE                                                          
!                                                                       
!-----  DAVIDON'S PROCEDURE  -----                                      
!     MINIMIZATION OF INNOVATION VARIANCE                               
!c      CALL MSDAV2(TL,SIGMA2,G,CR,Y,N,P0,IQ,IP,ISWRO,IPRNT)              
!x      CALL MSDAV2(TL2,SIGMA2,G2,CR,Y,N,P02,IQ,IP,ISWRO,HS,IPRNT,LU)
call msdav2(tl2,sigma2,g2,cr,y,n,p02,iq,ip,iswro,hs)
!                                                                       
!------------------------------                                         
!x      IF(IPRNT.EQ.0) GO TO 3930                                         
!c      WRITE(6,1200) ISWRO                                               
!x      WRITE(LU,1200) ISWRO                                               
!x 3930 CONTINUE                                                          
!                                                                       
if(iswro.ge.ipq) go to 1201
do 1902 i=1,ipq
!c      IF(DABS(PP0(I)-P0(I)).GE.CST005) GO TO 1919                       
if(dabs(pp0(i)-p02(i)).ge.cst005) go to 1919
1902 continue
go to 1201
1919 continue
!                                                                       
!x      IF(IPRNT.EQ.0) GO TO 3950                                         
!c      WRITE(6,1926)                                                     
!x      WRITE(LU,1926)                                                     
!x 3950 CONTINUE                                                          
!                                                                       
go to 4890
!                                                                       
1201 continue
!------------------------------                                         
!                                                                       
!c      WRITE(6,650)                                                      
!c      WRITE(6,520) (P0(I),I=1,IPQ)                                      
!c      WRITE(6,660)                                                      
!c      WRITE(6,520) (G(I),I=1,IPQ)                                       
!                                                                       
!-----  ALPH(I) PRINT OUT  -----                                        
!c      WRITE(6,6100)                                                     
icond=0
if(iq.le.0) go to 6510
do 6500 i=1,iq
!c 6500 PP0(I)=P0(I)                                                      
!c      CALL ARCHCK(PP0,IQ,ICOND)                                         
!xx 6500 PP0(I)=P02(I)                                                      
pp0(i)=p02(i)
6500 continue
call archck(pp0,alphb,iq,icond)
!c      WRITE(6,6501)                                                     
!c      WRITE(6,520) (ALPH(I),I=1,IQ)                                     
6510 if(ip.le.0) go to 6800
do 6700 i=1,ip
ii=iq+i
!c 6700 PP0(I)=P0(II)                                                     
!c      CALL ARCHCK(PP0,IP,ICOND)                                         
!xx 6700 PP0(I)=P02(II)
pp0(i)=p02(ii)
6700 continue
call archck(pp0,alpha,ip,icond)
!c      WRITE(6,6502)                                                     
!c      WRITE(6,520) (ALPH(I),I=1,IP)                                     
6800 continue
!------------------------------                                         
!                                                                       
!c      WRITE(6,670) TL                                                   
!                                                                       
!                                                                       
return
!xx  450 FORMAT(//1H ,5('-'),2X,'INITIAL VALUES P0(I)',2X,5('-'))          
!xx  500 FORMAT(//1H ,5('-'),2X,'INITIAL GRADIENT G(I)',2X,5('-'))         
!xx  520 FORMAT(1H ,5D20.10,/(1H ,5D20.10))                                
!xx  550 FORMAT(//1H ,'INITIAL (-2)LOG LIKELIHOOD=',D20.10,///)            
!xx  650 FORMAT(////1H ,5('-'),2X,'FINAL VALUES P0(I)',2X,5('-'))          
!xx  660 FORMAT(//1H ,5('-'),2X,'FINAL GRADIENT G(I)',2X,5('-'))           
!xx  670 FORMAT(//1H ,'FINAL (-2)LOG LIKELIHOOD=',D20.10,///)              
!xx 1200 FORMAT(//1H ,'ISWRO=',I5)                                         
!xx 1926 FORMAT(1H ,'HESSIAN RESET')                                       
!xx 3000 FORMAT(1H ,'INVERSE OF HESSIAN')                                  
!xx 3110 FORMAT(1H ,I5,4X,10D12.5,/(1H ,8X,10D12.5))                       
!xx 3900 FORMAT(1H ,'CR(I)')                                               
!xx 3910 FORMAT(1H ,5D20.10)                                               
!xx 6100 FORMAT(//1H ,5('-'),2X,'FINAL ALPH(I)''S AT SUBROUTINE ARCHCK',2X,
!xx     A5('-'))                                                           
!xx 6501 FORMAT(/1H ,12X,'AR-PART')                                        
!xx 6502 FORMAT(/1H ,12X,'MA-PART')                                        
end
subroutine subpm(p,b,a,iq,ip,ir)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!-----------------------------------------------------------------------
!     THIS SUBROUTINE COMPUTES THE VARIANCE MATRIX OF A STATIONARY STATE
!     VECTOR BY THE PROCEDURE OF AKAIKE (RESEARCH MEMO. 139 INST. STATIS
!     MATH. OCTOBER, 1978).                                             
!                                                                       
!     INPUTS:                                                           
!          IQ:  AR-ORDER                                                
!          (B(I),I=1,IQ):  AR-COEFFICIENTS                              
!          IP:  MA-ORDER                                                
!          (A(I),I=1,IP):  MA-COEFFICIENTS                              
!          IR:  IR=MAX(IQ,IP+1)                                         
!                                                                       
!     OUTPUTS:                                                          
!          ((P(I,J),I=1,IR),J=1,IR):  VARIANCE MATRIX OF THE STATIONARY 
!                                     STATE VECTOR                      
!-----------------------------------------------------------------------
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      DIMENSION P(51,51),A(50),B(50),WS(51),R(51),DB(1300)              
!xx      DIMENSION P(IR,IR),A(IR),B(IR),WS(IR),R(IR+1),DB(IQ*2)
integer iq, ip, ir
real(dp) p(ir,ir), b(ir), a(ir)
! local
integer i, ih, ihp1, ij1, ijk, im1, imik, imj, impr, ip1, ip2,&
&ipk, ipmi, ipp1, iqp1, irm1, irp1, it, itmk, itp1, itpk,&
&j, ji, k, l
real(dp) ws(ir), r(ir+1), db(iq*2), cst0, cst1, sum,&
&ckb, cki
data cst0,cst1/0.0d-00,1.0d-00/
!                                                                       
!                                                                       
!                                                                       
!c      DO 6100 I=1,1300                                                  
do 6100 i=1,iq*2
db(i)=cst0
6100 continue
!c      DO 6200 I=1,51                                                    
do 6200 i=1,ir
r(i)=cst0
6200 continue
!                                                                       
!-----  IMPULSE RESPONSE COMPUTATION  -----                             
ws(1)=cst1
ipp1=ip+1
iqp1=iq+1
irm1=ir-1
if(irm1.le.0) go to 400
do 100 i=2,ir
im1=i-1
impr=min0(iq,im1)
sum=cst0
if(impr.eq.0) go to 201
do 200 j=1,impr
imj=i-j
sum=sum-b(j)*ws(imj)
200 continue
201 if(i.le.ipp1) sum=sum+a(im1)
ws(i)=sum
100 continue
400 continue
!                                                                       
!-----  PREPARATION OF CONSTANT VECTOR  -----                           
irp1=ir+1
r(irp1)=cst0
if(irm1.eq.0) go to 7200
r(ir)=a(irm1)
do 7010 i=1,irm1
im1=i-1
ipmi=ip-im1
sum=cst0
if(ipmi.le.0) go to 7021
do 7020 j=1,ipmi
ji=im1+j
sum=sum+a(ji)*ws(j+1)
7020 continue
7021 if(i.eq.1) sum=sum+cst1
if(i.gt.1) sum=sum+a(im1)
r(i)=sum
7010 continue
go to 7300
7200 continue
r(1)=cst1
7300 continue
!                                                                       
!                                                                       
if(iq.eq.0) go to 1600
!                                                                       
!-----  TRIANGULARIZATION OF THE COEFFICIENT MATRIX  -----              
iqp1=iq+1
do 2100 l=1,iq
!xx 2100 DB(L)=B(L)
db(l)=b(l)
2100 continue
it=iq
i=iq
2200 continue
ip1=i+1
ip2=ip1+1
ih=ip2/2
ihp1=ih+1
ckb=db(it)
cki=cst1/(cst1-ckb**2)
do 2300 k=1,ih
ipk=ip2-k
r(k)=(r(k)-ckb*r(ipk))*cki
2300 continue
!                                                                       
if(i.le.2) go to 2401
do 2400 k=ihp1,i
ipk=ip2-k
r(k)=r(k)-ckb*r(ipk)
2400 continue
!                                                                       
2401 im1=i-1
!xx 2410 IF(IM1.EQ.0) GO TO 2600                                           
if(im1.eq.0) go to 2600
do 2500 k=1,im1
itpk=it+k
imik=it-i+k
itmk=it-k
db(itpk)=(db(imik)-ckb*db(itmk))*cki
2500 continue
!                                                                       
i=im1
it=it+i
go to 2200
2600 continue
!                                                                       
!-----  SOLVING THE LINEAR EQUATION  -----                              
if(iq.le.1) go to 3110
itp1=it+1
do 3100 i=2,iq
sum=r(i)
im1=i-1
do 3200 j=1,im1
itp1=itp1-1
sum=sum-db(itp1)*r(j)
3200 continue
r(i)=sum
3100 continue
3110 continue
do 3300 i=iqp1,irp1
sum=r(i)
do 3400 j=1,iq
imj=i-j
sum=sum-db(j)*r(imj)
3400 continue
r(i)=sum
3300 continue
!                                                                       
1600 continue
!                                                                       
!-----  P(I,J) COMPUTATION  -----                                       
do 9100 i=1,ir
do 9200 j=1,i
sum=cst0
if(j.eq.1) go to 9400
do 9300 k=1,j
ijk=k+i-j
sum=sum+ws(ijk)*ws(k)
9300 continue
9400 continue
ij1=i-j+1
p(i,j)=r(ij1)-sum
p(j,i)=p(i,j)
9200 continue
9100 continue
return
end
