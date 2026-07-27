! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine spgrh(y,n,lagh1,ifpl1,mode,period,cxx,cn,xmean,sd,aic,&
&parcor,pxx,ier)
  use timsac_kinds, only: dp
  implicit none
!
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      INTEGER H,H1                                                      
!xx      INTEGER PERIOD
!C      REAL*4 X,SXX,BMYA,BMYM                                            
!c      DIMENSION Y(1)                                                    
!c      DIMENSION X(4000),CXX(1001)                                       
!c      DIMENSION A(101),B(101)                                           
!c      DIMENSION SXX(501)                                                
!xx      DIMENSION Y(N)
!xx      DIMENSION CXX(LAGH1),CN(LAGH1)
!xx      DIMENSION A(IFPL1),B(IFPL1)
!xx      DIMENSION PXX(LAGH1)
!xx      DIMENSION SD(IFPL1),AIC(IFPL1),PARCOR(IFPL1-1)
integer n, lagh1, ifpl1, mode, period, ier
real(dp) y(n), cxx(lagh1), cn(lagh1), xmean, sd(ifpl1),&
&aic(ifpl1), parcor(ifpl1-1), pxx(lagh1)
! local
integer k, l
real(dp) a(ifpl1), b(ifpl1), sgme2, oaic
!
!c      H=60                                                              
!c      H1=H+1                                                            
!                                                                       
!c      DO 10 I=1,N                                                       
!c      X(I)=Y(I)                                                         
!c   10 CONTINUE                                                          
!                                                                       
!     AUTO COVARIANCE COMPUTATION.                                      
!c      LAGH=N-1                                                          
!c      LAGH=MIN0(LAGH,60)                                                
!c      LAGH1=LAGH+1                                                      
!c      CALL SAUTCO(X,CXX,N,LAGH1)
call autcorf(y,n,cxx,cn,lagh1,xmean)
!c      AN=N                                                              
!c      IFPL=3.0D-00*DSQRT(AN)                                            
!c      IFPL=MIN0(IFPL,50,LAGH)                                           
!c      IFPL=MIN0(30,N-1)                                                 
!c      IFPL1=IFPL+1                                                      
!c      CALL SICP2(CXX,IFPL1,N,A,L,SGME2,OAIC)                            
call sicp2(cxx,ifpl1,n,a,l,sgme2,oaic,sd,aic,parcor,ier)
if(mode .eq. 0) return
!c      WRITE(6,50)                                                       
k=0
!c 1031 CALL SNRASP(A,B,SXX,SGME2,L,K,H1,BMYA,BMYM)                       
!xx 1031 CALL SNRASP(A,B,PXX,SGME2,L,K,LAGH1,PERIOD)
call snrasp(a,b,pxx,sgme2,l,k,lagh1,period)
!
return
!xx   50 FORMAT(1H ,'*** HIGH ORDER AR-SPECTRUM AS AN APPROXIMATION TO PERI
!xx     *ODGRAM (ORDER IS FIXED AT 30 ) ***')                              
end
!c      SUBROUTINE SICP2(CYY,L1,N,COEF,MO,OSD,OAIC)                       
subroutine sicp2(cyy,l1,n,coef,mo,osd,oaic,sd1,aic1,parcor,ier)
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
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!xx      DIMENSION CYY(L1),COEF(L1)                                        
!c      DIMENSION A(101),B(101)                                           
!xx      DIMENSION A(L1),B(L1)
!xx      DIMENSION SD1(L1),AIC1(L1),PARCOR(L1-1)
integer l1, n, mo, ier
real(dp) cyy(l1), coef(l1), osd, oaic, sd1(l1), aic1(l1),&
&parcor(l1-1)
! local
integer i, ian, ian1, ian2, im, jj, jj0, jjl, l, lan1, lan2, lm,&
&m, mp1, nfc
real(dp) a(l1), b(l1), cst0, cst1, cst2, cst20, cst05,&
&cst01, sd, an, ran, scalh, aic, se, sdr, d, d2,&
&am, anfc

!c      REAL*4 AX,BL,STA,DASH,PLUS                                        
!c      REAL*4 FFFF                                                       
!c      REAL*4  F(41) / 41*1H  /, AMES(41) / 41*1H- /                     
!c      DATA AX,BL,STA,DASH,PLUS/1H!,1H ,1H*,1H-,1H+/                     
!-----
ier=0
!-----
cst0=0.0d-00
cst1=1.0d-00
cst2=2.0d-00
cst20=20.0d-00
cst05=0.05d-00
cst01=0.00001d-00
l=l1-1
sd=cyy(1)
an=n
oaic=an*dlog(sd)
osd=sd
mo=0
!     INITIAL CONDITION PRINT OUT                                       
!c  991 CONTINUE                                                          
!     WRITE(6,1100)                                                     
ran=cst1/dsqrt(an)
scalh=cst20
!xx      JJ0=SCALH+CST1                                                    
!xx      JJL=SCALH*CST2+CST1                                               
jj0=int(scalh+cst1)
jjl=int(scalh*cst2+cst1)
!c      AMES(1)=PLUS                                                      
!c      AMES(11)=PLUS                                                     
!c      AMES(JJ0)=PLUS                                                    
!c      AMES(JJ0+10)=PLUS                                                 
!c      AMES(JJL)=PLUS                                                    
!xx      IAN=SCALH*(RAN+CST05)                                             
ian=int(scalh*(ran+cst05))
ian1=ian+jj0
ian2=2*ian+jj0
lan1=-ian+jj0
lan2=-2*ian+jj0
!     WRITE(6,26100)                                                    
!c      WRITE(6,26101)                                                    
!c      WRITE(6,261)                                                      
!c      WRITE(6,26102)                                                    
!c      WRITE(6,262)                                                      
!c      WRITE(6,264) (AMES(J),J=1,JJL)                                    
!c      WRITE(6,859) MO,OSD,OAIC
sd1(mo+1)=osd
aic1(mo+1)=oaic
!-----
aic=oaic
!-----
!c      F(JJ0)=AX                                                         
!c      F(IAN1)=AX                                                        
!c      F(IAN2)=AX                                                        
!c      F(LAN1)=AX                                                        
!c      F(LAN2)=AX                                                        
!c      WRITE(6,861) (F(J),J=1,JJL)                                       
se=cyy(2)
!     ITERATION START                                                   
do 400 m=1,l
sdr=sd/cyy(1)
if(sdr.ge.cst01) go to 399
!c      WRITE(6,2600)                                                     
!-----
ier=2600
!-----
go to 402
399 mp1=m+1
d=se/sd
a(m)=d
d2=d*d
sd=(cst1-d2)*sd
am=m
aic=an*dlog(sd)+cst2*am
if(m.eq.1) go to 410
!     A(I) COMPUTATION                                                  
lm=m-1
do 420 i=1,lm
a(i)=a(i)-d*b(i)
420 continue
410 continue
do 421 i=1,m
im=mp1-i
!xx  421 B(I)=A(IM)
b(i)=a(im)
421 continue
!     M,SD,AIC  PRINT OUT                                               
if(a(m).lt.cst0) go to 300
!xx      NFC=SCALH*(A(M)+CST05)                                            
nfc=int(scalh*(a(m)+cst05))
go to 310
!xx  300 NFC=SCALH*(A(M)-CST05)                                            
300 nfc=int(scalh*(a(m)-cst05))
310 anfc=nfc
!xx      JJ=ANFC+SCALH+CST1                                                
jj=int(anfc+scalh+cst1)
!c      FFFF=F(JJ)                                                        
!c      F(JJ)=STA                                                         
!c      WRITE(6,860) M,SD,AIC,A(M),(F(J),J=1,JJL)                         
sd1(m+1)=sd
aic1(m+1)=aic
parcor(m)=a(m)
!c      F(JJ)=FFFF                                                        
!                                                                       
!xx  990 IF(OAIC.LT.AIC) GO TO 440                                         
if(oaic.lt.aic) go to 440
oaic=aic
osd=sd
mo=m
go to 440
!     DO 430 I=1,M                                                      
! 430 COEF(I)=-A(I)                                                     
440 if(m.eq.l) go to 400
se=cyy(m+2)
do 441 i=1,m
!xx  441 SE=SE-B(I)*CYY(I+1)                                               
se=se-b(i)*cyy(i+1)
441 continue
400 continue
402 continue
!c      WRITE(6,870) OAIC,MO                                              
oaic=aic
osd=sd
mo=l
do 5100 i=1,l
!xx 5100 COEF(I)=-A(I)
coef(i)=-a(i)
5100 continue
!                                                                       
!     MO, COEF(I) OUT PUT                                               
!     WRITE(6,1871)                                                     
!     WRITE(6,871)                                                      
!     CALL SUBVCP(COEF,MO)                                              
!c  699 F(JJ0)=BL                                                         
!c      F(IAN1)=BL                                                        
!c      F(IAN2)=BL                                                        
!c      F(LAN1)=BL                                                        
!c      F(LAN2)=BL                                                        
!c      AMES(JJ0)=DASH                                                    
!c      AMES(JJ0+10)=DASH                                                 
!c      AMES(JJL)=DASH                                                    
return
!xx26100 FORMAT(1H ,16X,'SD(M)',15X,'AIC(M)',13X,'A(M)')                   
!xx26101 FORMAT(1H ,////,4X,'M',11X,'INNOVATION',10X,'AIC(M)=    ',8X,     
!xx     A      'PARTIAL AUTO-')                                            
!xx  261 FORMAT(1H ,16X,'   VARIANCE',9X,'N*DLOG(SD(M))+2*M',              
!xx     A      '  CORRELATION    ',                                        
!xx     A      'PARTIAL CORRELATION (LINES SHOW +SD AND +2SD)')
!xx26102 FORMAT(1H+,102X,'_',7X,'_')                                       
!xx  262 FORMAT(1H ,69X,'-1',19X,'0',19X,'1')                              
!xx  264 FORMAT(1H ,70X,41A1)                                              
!xx  859 FORMAT(1H ,I5,2X,2D20.5)                                          
!xx  860 FORMAT(1H ,I5,2X,3D20.5,3X,41A1)                                  
!xx  861 FORMAT(1H+,70X,41A1)                                              
!xx  960 FORMAT(1H ,5X,'A(I)')                                             
!xx  870 FORMAT(/,1H ,'MINIMUM AIC(M)=',D12.5,2X,'ATTAINED AT M=',I5/)     
!xx  980 FORMAT(1H ,5X,'COEF(I)')                                          
!xx 1100 FORMAT(1H ,'AIC(M)=N*DLOG(SD)+2*M')                               
!xx  871 FORMAT(1H ,'AR-COEFFICIENTS')                                     
!xx 1871 FORMAT(1H ,'AR MODEL: Y(N)+AR(1)Y(N-1)+...+AR(M)Y(N-M)=X(N)')     
!xx 2600 FORMAT(1H ,'ACCURACY OF COMPUTATION LOST')                        
end
!c      SUBROUTINE SNRASP(A,B,SXX,SGME2,L,K,H1,BMYA,BMYM)                 
subroutine snrasp(a,b,pxx,sgme2,l,k,h1,ippp)
  use timsac_kinds, only: dp
  implicit none
!     THIS PROGRAM COMPUTES POWER SPECTRUM OF AN AR-MA PROCESS DEFINED B
!     X(N)=A(1)X(N-1)+...+A(L)X(N-L)+E(N)+B(1)E(N-1)+...+B(K)E(N-K),    
!     WHERE E(N) IS A WHITE NOISE WITH ZERO MEAN AND VARIANCE EQUAL TO  
!     SGME2.  OUTPUTS PXX(I) ARE GIVEN AT FREQUENCIES I/(2*H)           
!     I=0,1,...,H.                                                      
!     REQUIRED INPUTS ARE:                                              
!     L,K,H,SGME2,(A(I),I=1,L), AND (B(I),I=1,K).                       
!     0 IS ALLOWABLE AS L AND/OR K.                                     
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      INTEGER H,H1                                                      
!xx      INTEGER H1                                                      
!C      REAL*4 SXX,BMYA,BMYM,FY
!c      REAL*4 FY                                           
!c      INTEGER F,AX,BL,STA                                               
!c      DIMENSION F(61),FY(10)                                            
!c      DIMENSION A(101),B(101)                                           
!c      DIMENSION G(501),GR1(501),GI1(501),GR2(501),GI2(501)              
!c      DIMENSION PXX(501),SXX(501)                                       
!xx      DIMENSION A(L),B(K)                                           
!xx      DIMENSION G(L+K+1),GR1(H1),GI1(H1),GR2(H1),GI2(H1)
!xx      DIMENSION PXX(H1)
integer l, k, h1, ippp
real(dp) a(l), b(k), pxx(h1), sgme2
! local
integer i, i1, ippp0, k1, l1
real(dp) g(l+k+1), gr1(h1), gi1(h1), gr2(h1), gi2(h1),&
&cst0, cst1, t0
!c      COMMON /ILOGT/IDUMMY(2),IPUNC                                     
!c      COMMON /IDATA/IPPP                                                
!c      DATA AX/1H!/                                                      
!c      DATA BL,STA/1H ,1H*/                                              
ippp0=120/ippp
cst0=0.0d-00
cst1=1.0d-00
!xx  310 CONTINUE                                                          
if(l.le.0) go to 320
do 330 i=1,l
!xx  330 A(I)=-A(I)                                                        
a(i)=-a(i)
330 continue
320 l1=l+1
k1=k+1
g(1)=cst1
if(l.le.0) go to 400
do 10 i=1,l
i1=i+1
!xx   10 G(I1)=-A(I)
g(i1)=-a(i)
10 continue
!     COMMON SUBROUTINE CALL                                            
400 call fouger(g,l1,gr1,gi1,h1)
g(1)=cst1
if(k.le.0) go to 410
do 20 i=1,k
i1=i+1
!xx   20 G(I1)=B(I)                                                        
g(i1)=b(i)
20 continue
!     COMMON SUBROUTINE CALL                                            
410 call fouger(g,k1,gr2,gi2,h1)
do 30 i=1,h1
!xx   30 PXX(I)=(GR2(I)**2+GI2(I)**2)/(GR1(I)**2+GI1(I)**2)*SGME2
pxx(i)=(gr2(i)**2+gi2(i)**2)/(gr1(i)**2+gi1(i)**2)*sgme2
30 continue
!     WRITE(6,60)                                                       
!     WRITE(6,160)                                                      
!     WRITE(6,61) L,K,H                                                 
!     WRITE(6,164) SGME2                                                
!     WRITE(6,264)                                                      
!     WRITE(6,265)                                                      
if(l.le.0) go to 500
do 340 i=1,l
!xx  340 A(I)=-A(I)
a(i)=-a(i)
340 continue
!     WRITE(6,62)                                                       
!     COMMON SUBROUTINE CALL                                            
!     CALL SUBVCP(A,L)                                                  
500 if(k.le.0) go to 510
!c      WRITE(6,63)                                                       
!     COMMON SUBROUTINE CALL                                            
!c      CALL SUBVCP(B,K)                                                  
510 t0=cst0
!c      DO 520 I=1,H1                                                     
!c      T=PXX(I)                                                          
!c      IF(T.LT.T0) T=-T                                                  
!c  520 SXX(I)=DLOG10(T)                                                  
!     SEARCH  FOR THE MAXIMUM AND MINIMUM OF (SXX(I),I=0,LAGH)          
!c      CALL DSP3(SXX,H1,BMYA,BMYM)                                       
!c      WRITE(6,266)                                                      
!c  266 FORMAT(1H ,47X,'(DENSITY IN DB)')                                 
!c      FY(1)=BMYM*10.0                                                   
!c      DO 2031 I=2,7                                                     
!c 2031 FY(I)=FY(I-1)+10.0                                                
!c      WRITE(6,2033) (FY(I),I=1,7)                                       
!c 2033 FORMAT(1H ,25X,6(F5.1,5X),F5.1)                                   
!c      WRITE(6,67)                                                       
!c   67 FORMAT(1H ,4X,'I',1X,'RATIONAL SPECTRUM')                         
!c      WRITE(6,2035)                                                     
!c 2035 FORMAT(1H+,28X,6('+',9('-')),'+')                                 
!c      DO 2050 J=1,61                                                    
!c 2050 F(J)=BL                                                           
!c      DO 2040 I=1,H1                                                    
!c      IT=(SXX(I)-BMYM)*10.0+1.0                                         
!c      IM1=I-1                                                           
!c      F(1)=AX                                                           
!c      F(IT)=STA                                                         
!c      WRITE(6,2060) IM1,PXX(I),(F(J),J=1,61)                            
!c      MODP=MOD(IM1,IPPP0)                                               
!c      IF(MODP .EQ. 0) WRITE(6,1600)                                     
!c 1600 FORMAT(1H+,90X,'----X')                                           
!c      IF(IM1 .EQ. 5) WRITE(6,1601)                                      
!c 1601 FORMAT(1H+,98X,'X: IF PEAKS(TROUGHS) APPEAR')                     
!c      IF(IM1 .EQ. 6) WRITE(6,1602)                                      
!c 1602 FORMAT(1H+,101X,'AT THESE FRIQUENCIES,')                          
!c      IF(IM1 .EQ. 7) WRITE(6,1603)                                      
!c 1603 FORMAT(1H+,101X,'TRY LOWER(HIGHER) VALUES')                       
!c      IF(IM1 .EQ. 8) WRITE(6,1604)                                      
!c 1604 FORMAT(1H+,101X,'OF RIGID AND WATCH ABIC')                        
!c      IF(IM1 .EQ. 42 .AND. IPPP .EQ. 12) WRITE(6,1605)                  
!c 1605 FORMAT(1H+,90X,'----- IF A PEAK APPEARS HERE')                    
!c      IF(IM1 .EQ. 43 .AND. IPPP .EQ. 12) WRITE(6,1606)                  
!c 1606 FORMAT(1H+,101X,'TRY TRADING-DAY ADJUSTMENT')                     
!c      F(IT)=BL                                                          
!c 2040 CONTINUE                                                          
!c 2060 FORMAT(1H ,I5,2X,D14.5,7X,61A1)                                   
!c      IF(IPUNC .EQ. 0) RETURN                                           
!                                                                       
!c      WRITE(7,700) H1                                                   
!c      WRITE(7,701)                                                      
!c      WRITE(7,702) (SXX(I),I=1,H1)                                      
!c  700 FORMAT(I5)                                                        
!c  701 FORMAT('(6D12.5)')                                                
!c  702 FORMAT(6D12.5)                                                    
return
!xx   60 FORMAT(//1H ,'RATIONAL SPECTRUM')                                 
!xx   61 FORMAT(1H ,2HL=,I5,2X,2HK=,I5,2X,2HH=,I5)                         
!xx   62 FORMAT(1H ,5X,'AR(I)')                                            
!xx   63 FORMAT(1H ,5X,'MA(I)')                                            
!xx  160 FORMAT(1H ,17HINITIAL CONDITION)                                  
!xx  164 FORMAT(1H ,6HSGME2=,D12.5)                                        
!xx  264 FORMAT(1H ,'AR-MA MODEL:')                                        
!xx  265 FORMAT(1H ,'Y(N)+AR(1)Y(N-1)+...+AR(L)Y(N-L)=X(N)+MA(1)X(N-1)+...+
!xx     AMA(K)X(N-K)')                                                     
end
