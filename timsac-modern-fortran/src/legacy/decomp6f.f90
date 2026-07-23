! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine  decompf(data,n,ipar,trend,seasnl,ar,trad,noise,&
&para,imiss,omaxx,ier )
  use timsac_kinds, only: dp
  implicit none
!
!xx      PARAMETER (IOPT=1)
!xx      PARAMETER (NIP=9, NPA=26)
integer, parameter :: iopt=1
integer, parameter :: nip=9, npa=26
!xx      REAL*8    DATA(N), para(NPA), omaxx
!xx      INTEGER   IPAR(NIP)
!xx      REAL*8    TREND(N),SEASNL(N),AR(N),TRAD(N),NOISE(N)
!xx      INTEGER   PERIOD, SORDER
integer n, ipar(nip), imiss, ier
real(dp) data(n), trend(n), seasnl(n), ar(n), trad(n),&
&noise(n), para(npa), omaxx
! local
integer lm1
! common /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
! common /CCC/
integer isw, ismt, idif, log, mesh
!xx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH
common    /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
common    /ccc/     isw, ismt, idif, log, mesh
!
para(1:npa) = 0.0d0
!xx      CALL  SPARAM0( N,IPAR,NIP,para,NPA )
call  sparam0( ipar,nip )
lm1 = l+m+1
!
!xx      call decompff (DATA,N,IPAR,TREND,SEASNL,AR,TRAD,NOISE,
call decompff (data,n,trend,seasnl,ar,trad,noise,&
!x     *               para,iopt,imiss,omaxx,LM1)
&para,iopt,imiss,omaxx,lm1,ier)
!
return
end
!
!xx      SUBROUTINE  DECOMPFF(DATA,N,IPAR,TREND,SEASNL,AR,TRAD,NOISE,
subroutine  decompff(data,n,trend,seasnl,ar,trad,noise,&
!x     *                     para,iopt,imiss,omaxx,LM1 )
&para,iopt,imiss,omaxx,lm1,ier )
  use timsac_kinds, only: dp
  implicit none
!
!                                                                       
!  Bug fixed (97/10/17)
!
!  99/8/12
!  iopt < 0 -> use para. not optimized. 
!  input 
!  para(1):tau21, para(2):tau22, para(3):tau23,
!  para(4,5,6,...): ARCOEF
!
!  m4 = 0(no trade), 1(2-TDF model), 6(7-TDF model) 
!
!  99/9/2  if iopt < 0 , need ispan (=ipar(9))
! 

!  ...  TIME SERIES DECOMPOSITION (SEASONAL ADJUSTMENT) ...             
!                                                                       
!     THE BASIC MODEL:                                                  
!                                                                       
!          y(n) = T(n) + AR(n) + S(n) + TD(n) + R(n) + W(n)             
!                                                                       
!       where                                                           
!            T(n):       trend component                                
!            AR(n):      AR process                                     
!            S(n):       seasonal component                             
!            TD(n):      trading day factor                             
!            R(n):       any other explanetory variables                
!            W(n):       observational noise                            
!                                                                       
!     COMPONENT MODELS:                                                 
!                                                                       
!       Trend component                                                 
!            T(N) =  T(N-1) + V1(N)                             :M1 = 1 
!            T(N) = 2T(N-1) - T(N-2) + V1(N)                    :M1 = 2 
!            T(N) = 3T(N-1) -3T(N-2) + T(N-2) + V1(N)           :M1 = 3 
!                                                                       
!       AR componet:                                                    
!            AR(n) = a(1)AR(n-1) + ... + a(m2)AR(n-m2) + V2(n)          
!                                                                       
!       Seasonal component:                                             
!            S(N) =  -S(N-1) - ... - S(N-PERIOD+1) + V3(N)    :SORDER=1 
!            S(N) = -2S(N-1) - ... -PERIOB*S(N-PERIOD+1)      :SORDER=2 
!                            - ... - S(n-2PERIOD+2) + V3(n)             
!       Trading day effect:                                             
!            TD(n) = b(1)*TRADE(n,1) + ... + b(7)*TRADE(n,7)            
!                                                                       
!            TRADE(n,i):  number of i-th days of the week in n-th data  
!            b(1) + ... + b(7) = 0                                      
!                                                                       
!     REFERENCES                                                        
!                                                                       
!       G. KITAGAWA (1981), A Nonstationary Time Series Model and Its   
!            Fitting by a Recursive Filter, Journal of Time Series      
!            Analysis, Vol.2, 103-116.                                  
!                                                                       
!       W. GERSCH and G. KITAGAWA (1983), The prediction of time series 
!            with Trends and Seasonalities, Journal of Business and     
!            Economic Statistics, Vol.1, 253-264.                       
!                                                                       
!       G. KITAGAWA (1984), A smoothness priors-state space modeling of 
!            Time Series with Trend and Seasonality, Journal of American
!            Statistical Association, VOL.79, NO.386, 378-389.          
!                                                                       
!     STRUCTURE OF THE PROGRAM                                          
!                                                                       
!       <DECOMP>                                                        
!          |---<SPARAM>                                                 
!          |      |---<ID>                                              
!          |      +---<PARCOR>                                          
!          |---<REDATA>                                                 
!          |---<AREA>                                                   
!          |---<LOGTRF>                                                 
!          |---<TRADE>                                                  
!          |---<EPARAM>                                                 
!          |      |---<SETFGH>                                          
!          |      |---<OPTMIZ>                                          
!          |      |      +---<LINEAR>                                   
!          |      |             +---<FUNCND>                            
!          |      |                    +---<FUNCSA>                     
!          |      |                           |---<ARCOEF>              
!          |      |                           |---<SMOTH3>              
!          |      |                           |      |---<HUSHL7>       
!          |      |                           |      |---<HUSHL4>       
!          |      |                           |      |---<RECOEF>       
!          |      |                           |      +---<ID>           
!          |      |                           +---<STATE>               
!          |      +---<PPARA>                                          
!          |             +---<ARCOEF>                                   
!          |---<FUNCSA>                                                 
!          +---<PLOTDD>                                                 
!                 |---<MAXMIN>                                          
!                 |---<XYAXIS>                                          
!                 +---<PLOTD>                                           
!                                                                       
!     THE FOLLOWING CONTROL PARAMETERS ARE PRESET AS DEFAULT OPTION     
!          M1     = 2    :trend order(0, 1, 2 or 3)                     
!          M2     = 0    :AR order (less than 11, try 2 first)          
!          PERIOD = 12   :number of seasons in one period               
!          SORDER = 1    :seasonal order (0, 1 or 2)                    
!          TRADE  = 0    :trading day adjustment (if TRADE = 1)         
!          MT     = 1    :original data input device, MT = 5:card reader
!          BSPAN  = 300  :maximum data length in filtering              
!          ISPAN  = 100  :number of data for backward filtering         
!          MISING                                                       
!          TAU2(I)       :system noise variances (i=1,2,3)              
!          PAC(I)        :PARCOR (i=1,...,M2)                           
!          IPR    = 7    :print out control                             
!          IDIF   = 1    :numerical differencing (1 sided or 2 sided)   
!          LOG    = 0    :log transformation of data (if LOG = 1)       
!          YEAR          :the first year of the data                    
!          MESH   = 1    :draw mesh on the figure (if MESH > 0)         
!                                                                       
!     These options can be changed by using NAMELIST 'PARAM'            
!                                                                       
!          EXAMPLE:                                                     
!             &PARAM M2=2,LOG=1,&END                                    
!                                                                       
!                                                                       
!     -----  WRITTEN BY GENSHIRO KITAGAWA  ----END S                    
!                                                                       

!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!
!c      DIMENSION  DATA(N) ,IPAR(11) 
!c      real*4     title(20)
!c      REAL*8     A(40), YMEAN ,para(26)
!xx      PARAMETER (NIP=9, NPA=26)
integer, parameter :: nip=9, npa=26
!xx      DIMENSION  DATA(N) ,IPAR(NIP)
!xx      DIMENSION  DATA(N)
!xx      REAL*8     A(L+M2), YMEAN ,para(NPA)
!xx      REAL*8     TREND(N),SEASNL(N),AR(N),TRAD(N),NOISE(N)
!xx      INTEGER    IMIS(N), PERIOD, SORDER
!c      COMMON     /COMSM1/  WORK(300000)
!xx      DIMENSION  Z(N), E(L,LM1,N), TDAY(N,7)
integer n, iopt, imiss, lm1, ier
real(dp) data(n), trend(n), seasnl(n), ar(n), trad(n),&
&noise(n), para(npa), omaxx
! local
integer i, ifg, lll, imis(n)
! common /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
! common /CCC/
integer isw, ismt, idif, log, mesh
real(dp) a(l+m2), ymean, z(n), e(l,lm1,n), tday(n,7), ff
!c      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,
!c     *               NS, NI, MISING, IOUT, LL, NN, NYEAR,nmonth, 
!c     *                     NPREDS, NPREDE, IPRED 
!c      COMMON     /COMSM4/  NP1, NP2, NP3, NP4, NP5, NP6, NP7, NP8       
!      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG                    
!xx      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
!c      common     /cccout/  IMIS(3000)
common     /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
common     /ccc/     isw, ismt, idif, log, mesh
!                                                                       
!           ...  set control parameters  ...                            
!
!c      NN = N
!c      CALL  SPARAM( MT,A,IPAR,para,iopt )
lll = l+m2
!xx      CALL  SPARAM( N,A,LLL,IPAR,NIP,para,NPA,iopt )
call  sparam( a,lll,para,npa,iopt )
!      outmax = omaxx                      
do 123 i=1,n
imis(i) = 0
!c      if(data(i) .gt. omaxx)  imis(i) = 1
if (imiss .gt. 0) then
if(data(i) .gt. omaxx)  imis(i) = 1
else if (imiss .lt. 0) then
if(data(i) .lt. omaxx)  imis(i) = 1
end if
123 continue
!                                                                       
!           ...  read original data  ...                                
!                                                                       
!c      CALL  REDATA( MT,DATA,ISW,IPR,WORK,N,TITLE,YMEAN )                     
!xx      CALL  REDATAD( DATA,ISW,IPR,Z,N,YMEAN )
call  redatad( data,isw,z,n,ymean )
!                                                                       
!           ...  allocation of working area  ...                        
!                                                                       
!c      CALL  AREA( 300000,N )
!c      if(N .le. 0) then
!c      para(3) = -1.0
!c      return
!c      end if                                           
!                                                                       
!           ...  log transformation  ...                                
!                                                                       
!      IF(LOG .EQ. 1)  CALL  LOGTRF( WORK,N )                            
!c      CALL  LOGTRF( WORK,IMIS,N,LOG )
!x      CALL  LOGTRF( Z,IMIS,N,LOG )
call  logtrf( z,imis,n,log,ier )
if(ier .ne. 0) return
!                                                                       
!           ...  prepare calender for trading day adjustment  ...       
!
if( m4 .ne. 0)  then
!c      if( isea .eq. 12 )  CALL  TRADE( NYEAR,nmonth,N,WORK(NP2) )
!c      if( isea .eq. 4 )  CALL  TRADE2( NYEAR,nmonth,N,WORK(NP2) )
if( period .eq. 12 )  call  trade( nyear,nmonth,n,tday )
if( period .eq. 4 )  call  trade2( nyear,nmonth,n,tday )
endif
!                                                                       
!           ...  estimation of parameters  ...                          
!                                                                           

!c      CALL  EPARAM( A,TITLE ,iopt)
call  eparam( z,e,tday,imis,n,a,iopt)
!                                                                       
!      write(6,10)
! 10   format( '*****************UJHHBLB*************')
!           ...  smoothing with estimated parameters  ...               
!                                                                       
ismt = 1
!c      LLL = L+M2
!c      CALL  FUNCSA( 1,A,FF,IFG )                                        
call  funcsa( z,e,tday,imis,n,lm1,lll,a,ff,ifg)
!                                                                       
!           ...  plot estimated components  ...                         
!                                   

!c      call trpar( A,para)
call trpar( a,lll,para,npa)

!      para(2)= FF
!c      CALL PLOTDD(WORK(NP1),L+M+1,TITLE,A,WORK(NP2),WORK(NP3),WORK(NP7),                                   
!c     *           TREND,SEASNL,AR,TRAD,NOISE)
call plotdd(n,z,e,lm1,tday,trend,seasnl,ar,trad,noise)
!     CALL PRINT( WORK(NP1),L+M+1,WORK(NP7) )                           
!                                                                       
!	do 234 i=1,n
! 234	ar(i) = imis(i)
!
!       ar(1) = omaxx

return
end
!c      SUBROUTINE  ARCOEF( PAC,K,AR )                                    
subroutine  arcoefd( pac,k,ar )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  TRANSFORMATION FROM PARCOR TO AR COEFFICIENTS  ...              
!                                                                       
!       INPUTS:                                                         
!         PAC:  VECTOR OF PARTIAL AUTOCORRELATIONS                      
!         K:    ORDER OF THE MODEL                                      
!                                                                       
!       OUTPUTS:                                                        
!         AR:   VECTOR OF AR-COEFFICIENTS                               
!                                                                       
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  AR(K) , PAC(K) , W(100)
!xx      DIMENSION  AR(K) , PAC(K) , W(K)
integer k
real(dp) pac(k), ar(k)
! local
integer ii, im1, j, jj
real(dp) w(k)
!                                                                       
do  30     ii=1,k
ar(ii) = pac(ii)
w(ii)  = pac(ii)
im1 = ii - 1
if( im1 .le. 0 )     go to 30
do  10     j=1,im1
jj = ii - j
!xx   10 AR(J) = W(J) - PAC(II)*W(JJ)                                      
ar(j) = w(j) - pac(ii)*w(jj)
10 continue
if( ii .eq. k )     go to 40
do  20     j=1,im1
!   20 W(J) = PAC(J)  # modified  97/10/17                            
!xx   20 W(J) = ar(J)                                                     
w(j) = ar(j)
20 continue
30 continue
40 continue
return
end
!c      SUBROUTINE  EPARAM( A,TITLE ,iopt)                                     
subroutine  eparam( z,e,tday,imis,n,a,iopt)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  Estimation of parameters  ...                                   
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      DIMENSION  A(40), AI(20)                                          
!xx      DIMENSION  A(L+M2), AI(L+M2)
!xx      DIMENSION  Z(N), E(L,L+M+1,N), TDAY(N,7), IMIS(N)
!xx      INTEGER    PERIOD, SORDER
integer n, iopt, imis(n)
! common /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
real(dp) z(n), e(l,l+m+1,n), tday(n,7), a(l+m2)
! local
integer i, lm
real(dp) ai(l+m2)
! common /CCC/
integer isw, ismt, idif, log, mesh
! common /COMSM3/
real(dp) f1, f2, f3, a1, a2, a3, di, ui, tdf
! common /CMFUNC/
real(dp) djacob, fc, sig2, aic, fi, sig2i, aici, gi, gc
!c      REAL*4     CPUS, CPUE, TITLE(20), TIME(3)
!c      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,        
!c     *                  NS, NI, MISING, IOUT,LL,N,NYEAR,nmonth,NPREDS,      
!c     *                     NPREDE, NPRED
!cx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
!c      COMMON    /CMFUNC/  DJACOB,F,SIG2,AIC,FI,SIG2I,AICI,GI(20),G(20) 
!      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG                   
!cx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
!xx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH
common  /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
common  /comsm3/  f1(10), f2(10), f3(300), a1(10), a2(10),&
&a3(300) ,di, ui(3), tdf(7)
common  /cmfunc/  djacob,fc,sig2,aic,fi,sig2i,aici,gi(200),gc(200)
common  /ccc/     isw, ismt, idif, log, mesh
external   funcsa
!                                                                       
lm = l + m2
ismt = 0
!      CALL  CLOCK( CPUS,3 )                                             
do 10 i=1,lm
!xx   10 AI(I) = A(I)                                                      
ai(i) = a(i)
10 continue
!                                                                       
call  setfgh
if(iopt .ge. 0) then
!c      CALL  OPTMIZ( FUNCSA,A,L+M2 )
call  optmiz( funcsa,z,e,tday,imis,n,a,lm,l,l+m+1 )
end if
!                                                                       
!      CALL  CLOCK( CPUE,4 )                                             
!      CALL  DATE( DAY )                                                 
!      CALL  CLOCK( TIME,1 )                                             
!                                                                       
!      WRITE(6,650)                                                      
!      WRITE(6,600)  (TITLE(I),I=1,10), N                                
!      IF(LOG .EQ. 1)  WRITE(6,660)  DJACOB                              
!      WRITE(6,610)  M1, M2, M3, M4, M5                                  
!      WRITE(6,620)  DAY, (TIME(I),I=1,2)                                
!                                                                       
!      WRITE(6,630)                                                      
!c      CALL  PPARA( AI,GI,L,M2,FI,AICI,SIG2I )                          
!                                                                       
!      WRITE(6,640)                                                      
!c      CALL  PPARA( A,G,L,M2,F,AIC,SIG2 )                               
!c      CALL  PPARA( A,GC,L,M2,FC,AIC,SIG2 )                               
!c      IF( M4 .EQ. 0 )  GO TO 30                                         
!      WRITE(6,680)                                                      
!      DO 20 I=1,7                                                       
!   20 WRITE(6,690)  I, TDF(I)                                           
!c   30 CONTINUE                                                          
!      CPU = CPUE - CPUS                                                 
!      WRITE(6,670)   CPU                                                
!      WRITE(6,650)                                                      
!                                                                       
return
!xx  600 FORMAT( ////10X,'---  DATA  ---',/,10X,10A4,5X,'N =',I4 )
!xx  610 FORMAT( /10X,'---  MODEL  ---',/,10X,'M1 =',I2,5X,'M2 =',I2,5X,   
!xx     *        'M3 =',I3,5X,'M4 =',I2,5X,'M5 =',I2 )                     
!xx  620 FORMAT( /,10X,'---  PROGRAM  DECOMP  ---',/,10X,'DATE: ',A8,5X,   
!xx     *        'TIME: ',3A4 )                                            
!xx  630 FORMAT( ///15X,'<<<  INITIAL ESTIMATES  >>>' )                    
!xx  640 FORMAT( ///15X,'<<<  FINAL ESTIMATES  >>>' )                      
!xx  660 FORMAT( 10X,'***  LOG TRANSFORMED  ***   LOG-JACOBIAN =',         
!xx     *        F12.4 )                                                   
!xx  650 FORMAT( 1H1 )                                                     
!xx  670 FORMAT( //,10X,'CPU TIME =',F10.2 )                               
!xx  680 FORMAT( //,T20,'TRADING DAY',/,T15,'I',T22,'FACTOR' )             
!xx  690 FORMAT( 14X,I1,D16.8 )                                            
end
!c      SUBROUTINE  FUNCND( FUNCT,M,A,F,G,IFG )                           
subroutine  funcnd( funct,z,e,tday,imis,n,m,a,f,g,ifg,l,lm1 )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  FUNCTION EVALUATION AND NUMERICAL DIFFERENCING  ...             
!                                                                       
!xx      IMPLICIT   REAL*8( A-H,O-Z )
!xx      DIMENSION  Z(N), E(L,LM1,N), TDAY(N,7), IMIS(N)
!c      DIMENSION  A(M) , G(M) , B(20)
!xx      DIMENSION  A(M) , G(M) , B(M)
integer n, m, ifg, l, lm1, imis(n)
real(dp) z(n), e(l,lm1,n), tday(n,7), a(m), f, g(m)
! local
integer i, icnt, ii
real(dp) b(m), const, fb, ff
! common /CCC/
integer isw, ismt, idif, log, mesh
! common /CMFUNC/
real(dp) djacob, fc, sig2, aic, fi, sig2i, aici, gi, gc
!      COMMON     / CCC /  ISW , IPR, ISMT, IDIF, log
!xx      COMMON     /CCC/    ISW, IPR, ISMT, IDIF, LOG, MESH              
!cx      COMMON     /CMFUNC/ DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
common   /ccc/    isw, ismt, idif, log, mesh
common   /cmfunc/ djacob,fc,sig2,aic,fi,sig2i,aici,gi(200),gc(200)
!
data       icnt /0/
const = 0.0001d0
!
!c      CALL  FUNCT( M,A,F,IFG )                                          
call  funct( z,e,tday,imis,n,lm1,m,a,f,ifg )
fb = f
if( isw .ge. 1 )   return
!                                                                       
!      WRITE( 6,600 )   (A(I),I=1,M)                                     
do 10  i=1,m
!xx   10 B(I) = A(I)                                                       
b(i) = a(i)
10 continue
!                                                                       
do 30  ii=1,m
b(ii) = a(ii) + const
!c      CALL  FUNCT( M,B,FF,IFG )                                         
call  funct( z,e,tday,imis,n,lm1,m,b,ff,ifg )
if( idif .eq. 1 )  go to 20
b(ii) = a(ii) - const
!c      CALL  FUNCT( M,B,FB,IFG )
call  funct( z,e,tday,imis,n,lm1,m,b,fb,ifg )
20 g(ii) = (ff-fb)/(const*idif)
if( g(ii) .gt. 1.0d20 )  g(ii) = (f-fb)/const
if( g(ii) .lt.-1.0d20 )  g(ii) = (ff-f)/const
if( fb.gt.f .and. ff.gt.f )  g(ii) = 0.0d0
!xx   30 B(II) = A(II)                                                     
b(ii) = a(ii)
30 continue
!                                                                       
!      WRITE( 6,610 )   (G(I),I=1,M)                                     
do 40 i=1,m
!xx   40 GC(I) = G(I)                                                      
gc(i) = g(i)
40 continue
icnt = icnt + 1
if(icnt .gt. 1)  return
!                                                                       
aici  = aic
sig2i = sig2
fi    = fc
do 50 i=1,m
!xx   50 GI(I) = G(I)                                                      
gi(i) = g(i)
50 continue
return
!xx  600 FORMAT( 3X,'---  PARAMETER  ---',(/,3X,5D13.5) )                  
!xx  610 FORMAT( 3X,'---  GRADIENT  ---',(/,3X,5D13.5) )                   
end
!c      SUBROUTINE  FUNCSA( KK,A,FF,IFG )                                 
subroutine  funcsa( z,e,tday,imis,n,lm1,kk,a,ff,ifg )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  Initial setting, filtering and smoothing  ...                   
!                                                                       
!xx      IMPLICIT REAL*8( A-H,O-Z )                                        

!c      DIMENSION  Y(40), S(40,40), R(40,40), A(KK)
!c      REAL*8     ZZ(3000)
!c      dimension   imisr(3000)
!xx      DIMENSION  Z(N), E(L,LM1,N), TDAY(N,7), IMIS(N), A(KK)
!
!xx      DIMENSION  IMISR(N), ZZ(N), Y(M), R(LM1+1,LM1), S(M+1,M+1)
!xx      INTEGER    PERIOD, SORDER
integer n, lm1, kk, ifg, imis(n)
! common /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
real(dp) z(n), e(l,lm1,n), tday(n,7), a(kk), ff
! local
integer i, j, m12, mj, ni, imisr(n)
real(dp) zz(n), y(m), r(lm1+1,lm1), s(m+1,m+1), tau2,&
&sum
! common /COMSM3/
real(dp) f1, f2, f3, a1, a2, a3, di, ui, tdf
! common /CCC/
integer isw, ismt, idif, log, mesh
!c      COMMON    /COMSM1/  WORK(300000)
!c      COMMON    /COMSM2/  K1, K2, K3, K4, K5, M, L, ISEA, KSEA,
!c     *              NS, NI, MISING, IOUT, LL, N, NYEAR, nmonth,NPREDS,  
!c     *                    NPREDE, NPRED
!cx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
!c      COMMON    /COMSM4/  NP1, NP2, NP3, NP4, NP5, NP6, NP7, NP8
!      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG       
!xx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
!c      common    /cccout/  IMIS(3000)      
common  /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
common  /comsm3/  f1(10), f2(10), f3(300), a1(10), a2(10),&
&a3(300), di, ui(3), tdf(7)
common  /ccc/     isw, ismt, idif, log, mesh
!                                                                       
y(1:m) = 0.0d0
ifg = 0
!c      K = K1 + K2 + K3 + K4 + K5                                        
!c      K12 = K1 + K2
!c      LK1 = L + K + 1                                    
m12 = m1 + m2
!c      LM1 = L + M + 1                                                   
!c      MJ = 40                                                           
mj=m+1


!                                                                       
!  ...  Backward Filtering  ...                                         
!                                                                       
!c      if(ni .gt. n) ni=n
ni=n
do 10  i=1,ni
imisr(i) = imis(ni-i+1)
!c   10 ZZ(I) = WORK(NI-I+1)                                              
!xx   10 ZZ(I) = Z(NI-I+1)                                              
zz(i) = z(ni-i+1)
10 continue
do 20  i=1,l
tau2 = 0.5d0*(1.0d0 + dsin( a(i) ))  + 0.0001
if(tau2 .lt. 1.0d-20)  tau2 = 1.0d-20
!xx   20 UI(I) = 1.0D0/DSQRT( TAU2 )                                     
ui(i) = 1.0d0/dsqrt( tau2 )
20 continue
!c      IF( K2 .EQ. 0 )  GO TO 40
!c      DO 30  I=1,K2
if( m2 .eq. 0 )  go to 40
do 30  i=1,m2
!xx   30 Y(I) = 0.90D0*DSIN( A(L+I) )                                      
y(i) = 0.90d0*dsin( a(l+i) )
30 continue
!c      CALL  ARCOEF( Y,K2,A2 )                                           
call  arcoefd( y,m2,a2 )
!                                                                       
!xx   40 DO 50  I=1,MJ
!xx      DO 50  J=1,MJ                                                     
!c      R(I,J) = 0.0D0                                                    
!xx   50 S(I,J) = 0.0D0
!xx      DO 55 I = 1,LM1+1      
!xx         DO 55 J = 1,LM1
!xx   55 R(I,J) = 0.0D0
40 continue
s(1:mj,1:mj) = 0.0d0
r(1:lm1+1,1:lm1) = 0.0d0
di = 1.0d0
!c   70 IF( K2 .EQ. 0 )  GO TO 90                                         
!c      F2(1) = 1.0D0/A2(K2)                                              
!c      IF(K2.EQ.1)  GO TO 90                                             
!c      DO 80  I=2,K2                                                     
!c   80 F2(I) = -A2(I-1)/A2(K2)                                           
if( m2 .eq. 0 )  go to 90
f2(1) = 1.0d0/a2(m2)
if(m2.eq.1)  go to 90
do 80  i=2,m2
!xx   80 F2(I) = -A2(I-1)/A2(M2)                                           
f2(i) = -a2(i-1)/a2(m2)
80 continue
90 continue
!                                                                       
!c      CALL  SMOTH3( ZZ,R,S,WORK(NP1),WORK(NP2),WORK(NP3),               
!c     *              WORK(NP4),WORK(NP5),IMISR,NI,LK1,MJ,FF,0,0 )
!     *              WORK(NP4),WORK(NP5),WORK(NP6),NI,LK1,MJ,FF,0,0 )
call  smoth3( zz,r,s,e,tday,imisr,n,lm1,ff,0,0 )

!                                                                       
!  ...  Transformation of State Vector  ...                             
!                                                                       
!c      CALL  RECOEF( S,K,K,MJ,Y )                                        
!c      CALL  STATE( Y,A1,K1 )                                            
!c      CALL  STATE( Y(K1+1),A2,K2 )                                      
!c      CALL  STATE( Y(K12+1),A3,K3 )                                     
call  recoef( s,m,m,m+1,y )
!xxx      CALL  STATE( Y,A1,M1 )
!xxx      CALL  STATE( Y(M1+1),A2,M2 )
!xxx      CALL  STATE( Y(M12+1),A3,M3 )
if(m1 .ne. 0) call  state( y,a1,m1 )
if(m2 .ne. 0) call  state( y(m1+1),a2,m2 )
if(m3 .ne. 0) call  state( y(m12+1),a3,m3 )
!                                                                       
!c      DO 210  I=1,K                                                     
do 210  i=1,m
sum = 0.0d0
!c      DO 220  J=I,K                                                     
do 220  j=i,m
!xx  220 SUM = SUM + S(I,J)*Y(J)
sum = sum + s(i,j)*y(j)
220 continue
!c  210 S(I,K+1) = SUM
!xx  210 S(I,M+1) = SUM
s(i,m+1) = sum
210 continue
!                                                                       
!  ...  Forward Filtering and/or Smoothing  ...                         
!                                                                       
!c      CALL  SMOTH3( WORK,R,S,WORK(NP1),WORK(NP2),WORK(NP3),
!c     *              WORK(NP4),WORK(NP5),IMIS,N,LK1,MJ,FF,1,ISMT )  
!     *              WORK(NP4),WORK(NP5),WORK(NP6),N,LK1,MJ,FF,1,ISMT )       
call  smoth3( z,r,s,e,tday,imis,n,lm1,ff,1,ismt )
!                                                                       
ff = -ff
return
!                                                                       
end
subroutine  hushl4( x,mj1,n,k,m,isw )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!          HOUSEHOLDER TRANSFORMATION;   TYPE 4                         
!                                                                       
!xx      IMPLICIT  REAL*8( A-H,O-Z )                                       
!c      DIMENSION  X(MJ1,K), D(50)                                        
!xx      DIMENSION  X(MJ1,K), D(K)
integer mj1, n, k, m, isw
real(dp) x(mj1,k)
! local
integer ii, j
real(dp) d(k), tol, d0, d1, h, g, s
data     tol /1.0d-30/
!                                                                       
if( isw .eq. 1 )   go to 110
!xx      DO 100  II=M,K                                                    
do 101  ii=m,k
d0 = x(ii,ii)
d1 = x(n,ii)
h = d0**2 + d1**2
if( h .gt. tol )  go to 20
g = 0.0d0
go to  100
20 g = dsqrt( h )
if( d0 .ge. 0.0d0 )   g = -g
h = h - d0*g
d0 = d0 - g
d(ii) = d0
!                                                                       
if( ii .eq. k )  go to 100
do 60  j=ii+1,k
s = (d0*x(ii,j) + d1*x(n,j))/h
x(ii,j) = x(ii,j) - d0*s
x(n,j) = x(n,j) - d1*s
60 continue
100 x(ii,ii) = g
101 continue
return
!                                                                       
!                                                                       
110 continue
do 120  j=m,k-1
s = (d(j)*x(j,k)+x(n,j)*x(n,k))
h = -x(j,j)*d(j)
s = s/h
x(j,k) = x(j,k)-d(j)*s
!xx  120 X(N,K) = X(N,K)-X(N,J)*S
x(n,k) = x(n,k)-x(n,j)*s
120 continue
return
!                                                                       
end
subroutine  hushl7( x,d,mj1,k,m,ke )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     Householder Transformation,  TYPE 7                               
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!xx      DIMENSION  X(MJ1,K) , D(MJ1)                                      
integer mj1, k, m, ke
real(dp) x(mj1,k), d(mj1)
! local
integer i, ii, ii1, j, ne
real(dp) tol, h, g, f, s
!                                                                       
tol = 1.0d-30
!                                                                       
!xx      DO 100  II=1,KE                                                   
do 101  ii=1,ke
ne = max(m,ii) + 1
h = 0.0d00
do 10  i=ii,ne
d(i) = x(i,ii)
!xx   10       H = H + D(I)*D(I)
h = h + d(i)*d(i)
10 continue
if( h .gt. tol )  go to 20
g = 0.0d00
go to 100
20 g = dsqrt( h )
f = x(ii,ii)
if( f .ge. 0.0d00 )   g = -g
d(ii) = f - g
h = h - f*g
!                                                                       
!          FORM  (I - D*D'/H) * X, WHERE H = D'D/2                      
!                                                                       
ii1 = ii+1
do 30 i=ii1,ne
!xx   30 X(I,II) = 0.0D0
x(i,ii) = 0.0d0
30 continue
if( ii .eq. k )  go to 100
!                                                                       
do 60  j=ii1,k
s = 0.0d00
do 40  i=ii,ne
!xx   40       S = S + D(I)*X(I,J)                                         
s = s + d(i)*x(i,j)
40 continue
s = s/h
do 50  i=ii,ne
!xx   50      X(I,J) = X(I,J) - D(I)*S                                     
x(i,j) = x(i,j) - d(i)*s
50 continue
60 continue
100 x(ii,ii) = g
101 continue
!                                                                       
return
!                                                                       
end
integer function  id( k )
  use timsac_kinds, only: dp
  implicit none
!
!  ...  ID = 1:    IF K > 0
!       ID = 0:    OTHERWISE
!
integer k
id = 0
if( k .gt. 0 )  id = 1
return
end
!c      SUBROUTINE  LINEA1( FUNCT,X,H,RAM,EE,G,K,IG )                     
subroutine  linea1( funct,z,e,tday,imis,n,l,lm1,&
&x,h,ram,ee,g,k,ig )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  LINE SEARCH (WOLFE'S ALGORITHM)  ...                            
!                                                                       
!xx      IMPLICIT  REAL  *8 ( A-H,O-Z )                                    
!xx      INTEGER  RETURN,SUB                                               
!c      DIMENSION  X(K) , H(K) , X1(20)                                   
!c      DIMENSION  G(20)                                                  
!xx      DIMENSION  Z(N), E(L,LM1,N), TDAY(N,7), IMIS(N)
!xx      DIMENSION  X(K) , H(K) , G(K), X1(K)
integer n, l, lm1, k, ig, imis(n)
real(dp) z(n), e(l,lm1,n), tday(n,7), x(k), h(k), ram,&
&ee, g(k)
! local
integer i, ifg, return, sub
real(dp) a1, a2, a3, b1, b2, c1, c2, e1, e2, e3, f0, h1,&
&h2, const2, hnorm, ram1, ram2, ram3, sum, sum0,&
&x1(k)
! common /CCC/
integer isw, ismt, idif, log, mesh
!      COMMON     / CCC /  ISW , IPR                          
!      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG
!xx      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH
common     /ccc/     isw, ismt, idif, log, mesh
external   funct
c1 = 0.01d0
c2 = 0.5d0
!
!c
ifg = 0
!c
f0 = ee
sum0 = 0.0d0
do 15 i=1,k
!xx   15 SUM0 = SUM0 + G(I)*H(I)                                           
sum0 = sum0 + g(i)*h(i)
15 continue
isw = 1
ram = 0.5d0
const2 = 1.0d-60
hnorm = 0.d0
do 10  i=1,k
!xx   10 HNORM = HNORM + H(I)**2                                           
hnorm = hnorm + h(i)**2
10 continue
hnorm = dsqrt( hnorm )
!                                                                       
ram2 = ram
e1 =ee
ram1 = 0.d0
!                                                                       
if( ram2*hnorm .gt. 5.0d0 )  go to 48
do 20  i=1,k
!xx   20 X1(I) = X(I) + RAM2*H(I)                                          
x1(i) = x(i) + ram2*h(i)
20 continue
!c      CALL  FUNCND( FUNCT,K,X1,E2,G,IG )                                
call  funcnd( funct,z,e,tday,imis,n,k,x1,e2,g,ig,l,lm1 )
!      IF(IPR.GE.7)  WRITE(6,2)  RAM2,E2                                 
!                                                                       
if( ig .eq. 1 )  go to  50
if( e2 .gt. e1 )  go to 50
!CCCCCCCCCCCCCCCCCCCCCC                                                 
h1 = e2 - f0 - c1*ram2*sum0
sum = 0.0d0
do 25 i=1,k
!xx   25 SUM = SUM + G(I)*H(I)                                             
sum = sum + g(i)*h(i)
25 continue
h2 = sum - c2*sum0
if(h1.le.0.0d0 .and. h2.ge.0.0d0)  return
!CCCCCCCCCCCCCCCCCCCCCC                                                 
30 ram3 = ram2*4.d0
do 40  i=1,k
!xx   40 X1(I) = X(I) + RAM3*H(I)
x1(i) = x(i) + ram3*h(i)
40 continue
!c      CALL  FUNCND( FUNCT,K,X1,E3,G,IG )                                
call  funcnd( funct,z,e,tday,imis,n,k,x1,e3,g,ig,l,lm1 )
if( ig.eq.1 )  go to  500
!      IF( IPR.GE.7 )  WRITE(6,3)  RAM3,E3                               
if( e3 .gt. e2 )  go to 70
if(ram3.gt.1.0d10 .and. e3.lt.e1)  go to 45
if(ram3.gt.1.0d10 .and. e3.ge.e1)  go to 46
ram1 = ram2
ram2 = ram3
e1 = e2
e2 = e3
go to 30
!                                                                       
45 ram = ram3
ee = e3
return
!                                                                       
46 ram = 0.0d0
return
!                                                                       
48 e2 = 1.0d30
!                                                                       
50 ram3 = ram2
e3 = e2
ram2 = ram3*0.1d0
if( ram2*hnorm .lt. const2 )  go to  400
do 60  i=1,k
!xx   60 X1(I) = X(I) + RAM2*H(I)                                          
x1(i) = x(i) + ram2*h(i)
60 continue
!c      CALL  FUNCND( FUNCT,K,X1,E2,G,IG )                                
call  funcnd( funct,z,e,tday,imis,n,k,x1,e2,g,ig,l,lm1 )
!      IF(IPR.GE.7)  WRITE(6,4)  RAM2,E2                                 
if( e2.gt.e1 )  go to 50
!CCCCCCCCCCCCCCCCCCCCCC                                                 
h1 = e2 - f0 - c1*ram2*sum0
sum = 0.0d0
do 65 i=1,k
!xx   65 SUM = SUM + G(I)*H(I)                                             
sum = sum + g(i)*h(i)
65 continue
h2 = sum - c2*sum0
if(h1.gt.0.0d0 .or. h2.lt.0.0d0)  go to 70
ram = ram2
return
!CCCCCCCCCCCCCCCCCCCCCC                                                 
!                                                                       
!c   70 ASSIGN 80 TO RETURN                                               
70 return = 80
go to 200
!                                                                       
80 do 90  i=1,k
!xx   90 X1(I) = X(I) + RAM*H(I)                                           
x1(i) = x(i) + ram*h(i)
90 continue
!c      CALL  FUNCND( FUNCT,K,X1,EE,G,IG )                                
call  funcnd( funct,z,e,tday,imis,n,k,x1,ee,g,ig,l,lm1 )
!      IF(IPR.GE.7)  WRITE(6,5)  RAM,EE                                  
!                                                                       
ifg = 0
!c      ASSIGN  300 TO  SUB                                               
!c      ASSIGN 200 TO SUB                                                 
!c   95 ASSIGN 130 TO RETURN                                              
sub = 300
sub = 200
95 return = 130
if( ram .gt. ram2 )  go to 110
if( ee .ge. e2 )  go to 100
ram3 = ram2
ram2 = ram
e3 =e2
e2 =ee
!c      GO TO  SUB,( 200,300 )                                            
if( sub .eq. 200 ) go to 200
if( sub .eq. 300 ) go to 300
!                                                                       
100 ram1 = ram
e1 = ee
!c      GO TO  SUB,( 200,300 )                                            
if( sub .eq. 200 ) go to 200
if( sub .eq. 300 ) go to 300
!                                                                       
110 if( ee .le. e2 )  go to 120
ram3 = ram
e3 = ee
!c      GO TO  SUB,( 200,300 )                                            
if( sub .eq. 200 ) go to 200
if( sub .eq. 300 ) go to 300
!                                                                       
120 ram1 = ram2
ram2 = ram
e1 = e2
e2 = ee
!c      GO TO  SUB,( 200,300 )                                            
if( sub .eq. 200 ) go to 200
if( sub .eq. 300 ) go to 300
!                                                                       
130 do 140  i=1,k
!xx  140 X1(I) = X(I) + RAM*H(I)                                           
x1(i) = x(i) + ram*h(i)
140 continue
!c      CALL  FUNCND( FUNCT,K,X1,EE,G,IG )                                
call  funcnd( funct,z,e,tday,imis,n,k,x1,ee,g,ig,l,lm1 )
!      IF( IPR.GE.7 )  WRITE(6,6)  RAM,EE                                
!c      ASSIGN 200 TO SUB                                                 
sub = 200
ifg = ifg+1
!     IFG = 0                                                           
!                                                                       
h1 = ee - f0 - c1*ram*sum0
sum = 0.0d0
do 145 i=1,k
!xx  145 SUM = SUM + G(I)*H(I)                                             
sum = sum + g(i)*h(i)
145 continue
h2 = sum - c2*sum0
if( h1.le.0.0d0 .and. h2.le.0.0d0 )  go to 150
if( ifg .le. 20 )  go to 95
!                                                                       
150 if( e2 .lt. ee )  ram = ram2
return
!                                                                       
!      -------  INTERNAL SUBROUTINE SUB1  -------                       
200 if( ram3-ram2 .gt. 5.0d0*(ram2-ram1) )  go to 202
if( ram2-ram1 .gt. 5.0d0*(ram3-ram2) )  go to 204
a1 = (ram3-ram2)*e1
a2 = (ram1-ram3)*e2
a3 = (ram2-ram1)*e3
b2 = (a1+a2+a3)*2.d0
b1 = a1*(ram3+ram2) + a2*(ram1+ram3) + a3*(ram2+ram1)
if( b2 .eq. 0.d0 )  go to 210
ram = b1 /b2
!c      GO TO RETURN ,( 80,130 )                                          
if( return .eq. 80 ) go to 80
if( return .eq. 130 ) go to 130
202 ram = (4.0d0*ram2 + ram3)/5.0d0
!c      GO TO RETURN, (80,130)                                            
if( return .eq. 80 ) go to 80
if( return .eq. 130 ) go to 130
204 ram = (ram1 + 4.0d0*ram2)/5.0d0
!c      GO TO RETURN, (80,130)                                            
if( return .eq. 80 ) go to 80
if( return .eq. 130 ) go to 130
!                                                                       
210 ig = 1
ram = ram2
return
!                                                                       
!      -------  INTERNAL SUBROUTINE SUB2  -------                       
!                                                                       
300 if( ram3-ram2 .gt. ram2-ram1 )  go to 310
ram = (ram1+ram2)*0.5d0
!c      GO TO RETURN ,( 80,130 )                                          
if( return .eq. 80 ) go to 80
if( return .eq. 130 ) go to 130
!                                                                       
310 ram = (ram2+ram3)*0.5d0
!c      GO TO RETURN ,( 80,130 )                                          
if( return .eq. 80 ) go to 80
if( return .eq. 80 ) go to 130
! ---------------------------------
!                                                                       
400 ram = 0.d0
return
! ------------------------------------------------------------          
!                                                                       
500 ram = (ram2+ram3)*0.5d0
510 do 520  i=1,k
!xx  520 X1(I) = X(I) + RAM*H(I)                                           
x1(i) = x(i) + ram*h(i)
520 continue
!c      CALL  FUNCND( FUNCT,K,X1,E3,G,IG )                                
call  funcnd( funct,z,e,tday,imis,n,k,x1,e3,g,ig,l,lm1 )
!      IF( IPR.GE.7 )  WRITE(6,7)  RAM,E3                                
if( ig.eq.1 )  go to 540
if( e3.gt.e2 )  go to 530
ram1 = ram2
ram2 = ram
e1 = e2
e2 = e3
go to 500
!                                                                       
530 ram3 = ram
go to 70
!                                                                       
540 ram = (ram2+ram)*0.5d0
go to 510
!                                                                       
! ------------------------------------------------------------          
!xx    1 FORMAT( 5X ,'RAM =',D10.3, 5X,'E1 =',F10.3 )                      
!xx    2 FORMAT( 5X ,'RAM =',D10.3, 5X,'E2 =',F10.3 )                      
!xx    3 FORMAT( 5X ,'RAM =',D10.3, 5X,'E3 =',F10.3 )                      
!xx    4 FORMAT( 5X ,'RAM =',D10.3, 5X,'E4 =',F10.3 )                      
!xx    5 FORMAT( 5X ,'RAM =',D10.3, 5X,'E5 =',F10.3 )                      
!xx    6 FORMAT( 5X ,'RAM =',D10.3, 5X,'E6 =',F10.3 )                      
!xx    7 FORMAT( 5X ,'RAM =',D10.3, 5X,'E7 =',F10.3 )                      
end
!c      SUBROUTINE  LOGTRF( Z,N,ilog )                                         
!x      SUBROUTINE  LOGTRF( Z,IMIS,N,ilog )
subroutine  logtrf( z,imis,n,ilog,ier )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
! ... LOG TRANSFORMATION ...                                            
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!xx      REAL*8   Z(N)
!xx      DIMENSION   IMIS(N)                                               
integer n, ilog, ier, imis(n)
real(dp) z(n)
! local
integer i
! common /CMFUNC/
real(dp) djacob, fc, sig2, aic, fi, sig2i, aici, gi, gc
!      COMMON     /CMFUNC/  DJACOB                                       
!c      COMMON     /CMFUNC/  DJACOB,F,SIG2,AIC,FI,SIG2I,AICI,GI(20),G(20)
!c      common     /cccout/  IMIS(3000)
!cx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
common  /cmfunc/  djacob,fc,sig2,aic,fi,sig2i,aici,gi(200),gc(200)
!                                                                       
ier = 0
djacob = 0.0d0
if(ilog .eq. 0) return
ier = -1
do 10 i=1,n
!      DJACOB = DJACOB - DLOG10(Z(I))                                    
!x	if(IMIS(i) .ne. 1)  DJACOB = DJACOB - DLOG(Z(I)) 
!   10 Z(I) = DLOG10( Z(I) )   
!x	Z(I) = DLOG( Z(I) )                                             
if(imis(i) .ne. 1) then
if(z(i) .le. 0) return
djacob = djacob - dlog(z(i))
z(i) = dlog( z(i) )
end if
10 continue
ier = 0
!                                                                       
!      WRITE(6,600)  DJACOB                                              
!     WRITE(6,610) (Z(I),I=1,N)                                         
return
!                                                                       
!xx  600 FORMAT(1H0,'OPTION = 1, (LOG TRANSFORMATION),     LOG-JACOBIAN =',
!xx     1       D15.7)                                                     
!xx  610 FORMAT(1H ,10D13.6)                                               
end
!c      SUBROUTINE  OPTMIZ( FUNCT,X,N )                                   
subroutine  optmiz( funct,z,e,tday,imis,nn,x,n,l,lm1 )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  NUMERICAL OPTIMIZATION  ...                                     
!            LATEST REVISION:  JUNE 20, 1983                            
!                                                                       
!xx      IMPLICIT  REAL*8 (A-H,O-Z)
!c      DIMENSION  X(20) , DX(20) , G(20) , G0(20) , Y(20)                
!c      DIMENSION  H(20,20) , WRK(20) , S(20)
!xx      DIMENSION  Z(NN), E(L,LM1,NN), TDAY(NN,7), IMIS(NN), X(N)
!xx      DIMENSION  DX(N) ,G(N) ,G0(N) ,Y(N) ,H(N,N), WRK(N), S(N)
integer nn, n, l, lm1, imis(nn)
real(dp) z(nn), e(l,lm1,nn), tday(nn,7), x(n)
! local
integer i, icc, icount, ig, j
real(dp) dx(n), g(n), g0(n), y(n), h(n,n), wrk(n), s(n),&
&tau2, eps1, eps2, const1, xm, sum, s1, s2, stem,&
&ss, ds2, gtem, ed, ramda, xmb
! common /CCC/
integer isw, ismt, idif, log, mesh
! common /CMFUNC/
real(dp) djacob, fc, sig2, aic, fi, sig2i, aici, gi, gc
!      COMMON     / CCC /  ISW, IPR                                      
!      COMMON     /CMFUNC/  DJACOB, F, SD, AIC                           
!c      COMMON     /CMFUNC/  DJACOB,F,SIG2,AIC,FI,SIG2I,AICI,
!c     *                   GI(20),GDUM(20) 
!xx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
!cx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
common  /ccc/     isw, ismt, idif, log, mesh
common  /cmfunc/  djacob,fc,sig2,aic,fi,sig2i,aici,gi(200),gc(200)
external   funct
!c      DATA  TAU1 , TAU2  /  1.0D-1 , 1.0D-1  /             
data  tau2  /  1.0d-1  /
data  eps1 , eps2  / 1.0d-2 , 1.0d-2  /
const1 = 1.0d-50
isw = 0
!                                                                       
!c      CALL  FUNCND( FUNCT,N,X,XM,G,IG )
xm = 0
call  funcnd( funct,z,e,tday,imis,nn,n,x,xm,g,ig,l,lm1 )
!                                                                       
!          INITIAL ESTIMATE OF INVERSE OF HESSIAN                       
!                                                                       
icount=0

!xx 1000 DO 20  I=1,N                                                      
!xx      DO 10  J=1,N                                                      
!xx   10 H(I,J) = 0.0D0                                                    
!xx      S(I)   = 0.0D0                                                    
!xx      DX(I)  = 0.0D0                                                    
!xx   20 H(I,I) = 1.0D0
1000 continue
h(1:n,1:n) = 0.0d0
s(1:n)     = 0.0d0
dx(1:n)    = 0.0d0
do 20  i=1,n
h(i,i) = 1.0d0
20 continue
icc = 0

icount=icount+1
if(icount .gt. 4) then
!         WRITE(6,700)
!         WRITE(6,701)
return
endif
!                                                                       
!      IF( IPR .GE. 2 )   WRITE( 6,340 )   XM                            
!                                                                       
!  ...  QUASI NEWTON ALGORITHM  ...                                     
!                                                                       
2000 icc = icc + 1
if( icc .eq. 1 )   go to 120
!                                                                       
do 40  i=1,n
!xx   40 Y(I) = G(I) - G0(I)                                               
y(i) = g(i) - g0(i)
40 continue
do 60  i=1,n
sum = 0.0d0
do 50  j=1,n
!xx   50 SUM = SUM + Y(J)*H(I,J)                                           
sum = sum + y(j)*h(i,j)
50 continue
!xx   60 WRK(I) = SUM                                                      
wrk(i) = sum
60 continue
s1 = 0.0d0
s2 = 0.0d0
do 70  i=1,n
s1 = s1 + wrk(i)*y(i)
!xx   70 S2 = S2 + DX(I) *Y(I)                                             
s2 = s2 + dx(i) *y(i)
70 continue
if( s1.le.const1 .or. s2.le.const1 )  go to 900
!                                                                       
!  ...  BFGS FORMULA FOR UPDATING INVERSE OF HESSIAN MATRIX  ...        
!                                                                       
stem = s1 / s2 + 1.0d00
!xx      DO 110  I=1,N                                                     
do 111  i=1,n
do 110  j=i,n
h(i,j) = h(i,j)- (dx(i)*wrk(j)+wrk(i)*dx(j)-dx(i)*dx(j)*stem)/s2
!xx  110 H(J,I) = H(I,J)                                                   
h(j,i) = h(i,j)
110 continue
111 continue
!                                                                       
120 ss = 0.0d0
do 150  i=1,n
sum = 0.0d0
do 140  j=1,n
!xx  140 SUM = SUM + H(I,J)*G(J)                                           
sum = sum + h(i,j)*g(j)
140 continue
ss = ss + sum * sum
!xx  150 S(I) = -SUM                                                       
s(i) = -sum
150 continue
!                                                                       
s1 = 0.0d0
s2 = 0.0d0
do 170  i=1,n
s1 = s1 + s(i)*g(i)
!xx  170 S2 = S2 + G(I)**2                                                 
s2 = s2 + g(i)**2
170 continue
ds2 = dsqrt(s2)
gtem = dabs(s1) / ds2
!     IF( GTEM .LE. TAU1  .AND.  DS2 .LE. TAU2 )     GO TO  900         
if( s1 .ge. 0.0d0 )  go to  1000
ed = xm
!                                                                       
!  ...  LINE SEARCH  ...                                                
!                                                                       
!c      CALL  LINEA1( FUNCT,X,S,RAMDA,ED,G,N,IG )                         
call linea1(funct,z,e,tday,imis,nn,l,lm1,x,s,ramda,ed,g,n,ig)
!                                                                       
!      IF( IPR .GE. 2 )   WRITE( 6,330 )   RAMDA, ED                     
!                                                                       
s1 = 0.0d0
do 210  i=1,n
dx(i) = s(i)*ramda
s1 = s1 + dx(i)**2
g0(i) = g(i)
!xx  210 X(I) = X(I) + DX(I)                                               
x(i) = x(i) + dx(i)
210 continue
xmb  = xm
isw  = 0
!                                                                       
!c      CALL  FUNCND( FUNCT,N,X,XM,G,IG )
call  funcnd( funct,z,e,tday,imis,nn,n,x,xm,g,ig,l,lm1 )
!                                                                       
s2 = 0.d0
do 220  i=1,n
!xx  220 S2 = S2 + G(I)**2
s2 = s2 + g(i)**2
220 continue
if( dsqrt(s2) .lt. tau2 )   go to 900
if( xmb-xm .lt. eps1  .and.  dsqrt(s1) .lt. eps2 )  go to 900
if( xmb-xm .lt. 0.0001 .and. icc .gt. n ) go to 900
!      IF( ICC .GE. N*2 )  GO TO 1000                                    
go to 2000
!                                                                       
900 continue
s2 = 0.0d0
do 230 i=1,n
!xx  230 S2 = S2 + G(I)**2                                                 
s2 = s2 + g(i)**2
230 continue
if( dsqrt(s2) .gt. 1.0d0 )  go to 1000
!xx      IF( IPR .LE. 0 )   RETURN                                         
!      WRITE(6,620)                                                      
!      WRITE( 6,600 )                                                    
!      WRITE( 6,610 )     (X(I),I=1,N)                                   
!      WRITE( 6,601 )                                                    
!      WRITE( 6,610 )     (G(I),I=1,N)                                   
return
!xx  330 FORMAT( 5X ,'RAM =',D10.3,5X,'E0 =',F10.3 )                       
!xx  340 FORMAT( 25X,'EE =',F10.3 )                                        
!xx  600 FORMAT( 1H0,'--  PARAMETER  ---' )                                
!xx  601 FORMAT( 1H0,'--  GRADIENT  --' )                                  
!xx  610 FORMAT( 1H ,10D13.5 )                                             
!xx  620 FORMAT( 10X,'<<<  FINAL ESTIMATES  >>>' )                         
!xx 700  FORMAT( ' OPTIMIZATION DOES NOT SUCCEED ' )
!xx 701  FORMAT( ' CHECK THE MODEL ! ' )
end

!c      SUBROUTINE  PLOTDD( E,LM1,TITLE,A,TRADE,REG,Z ,COMP1,COMP2,COMP3,
subroutine  plotdd( n,z,e,lm1,trade,comp1,comp2,comp3,&
&comp4,comp5)
  use timsac_kinds, only: dp
  implicit none
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          

!  ...  PLOT ESTIMATED PARAMETRS, ORIGINAL DATA AND ESTIMATED PARAMETERS
!
!xx      INTEGER    PERIOD, SORDER
integer n, lm1
!ommon /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
real(dp) z(n), e(l,lm1,n), trade(n,7), comp1(n), comp2(n),&
&comp3(n), comp4(n), comp5(n)
! local
integer i, id, j, m12, m123, m1234
real(dp) reg(n,m5), sum, tmp
! common /CCC/
integer isw, ismt, idif, log, mesh
!      COMMON   /COMSM1/  nwww, WORK                                        
!c      COMMON     /COMSM1/  WORK(300000)
!c      COMMON     /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,                  
!c     *               NS,NI,MISING,IOUT,LL,N,NYEAR,nmonth,NPS,NPE,NPRED
!xx      COMMON     /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
common     /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
common     /ccc/     isw, ismt, idif, log, mesh
!c      DIMENSION  E(L,LM1,N), Z(N), TRADE(N,7), REG(N,M5)
!xx      DIMENSION  E(L,LM1,N), Z(N), TRADE(N,7), REG(N,M5)
!xx      REAL*8     COMP1(N),COMP2(N),COMP3(N),COMP4(N),COMP5(N)
!c      real*4     TITLE(20)
!c      REAL*8     A(40),COMP1(N),COMP2(N),COMP3(N) 
!c      REAL*8     COMP4(N),COMP5(N)


!      CALL  PLTCK                                                       
!      CALL  PLOTS( 8HDATA,0.0 )                                         
!      CALL  NEWPEN( 2 )                                                 
!                                                                       
!      CALL  MAXMIN( WORK,N,Y0,Y1,DY )                                   
!                                                                       
!      XX = 1.5                                                          
!      YY = 15.0                                                         
!      X0 = 0                                                            
!      X1 = N                                                            
!      DX = 12.0                                                         
!      IF( N .GT. 200 )  DX = 24.0                                       
!      IF( N .GT. 400 )  DX = 48.0                                       
!      NS = 1                                                            
!      WX = N*0.1                                                        
!  150 IF(WX .LE. 13.0)  GO TO 160                                       
!      WX = WX *0.8                                                      
!      GO TO 150                                                         
!  160 CONTINUE                                                          
!      WY = 6.0                                                          
!      WC = 0.22                                                         
m12 = m1 + m2
m123= m1 + m2 + m3
m1234=m1 + m2 + m3 + m4
!                                                                       
!      NNY = (Y1-Y0)/DY                                                  
!      IF( NNY .GE. 4 )  DY1 = DY                                        
!      IF( NNY .EQ. 2 .OR. NNY .EQ. 3 )  DY1 = DY*0.5                    
!      IF( NNY .EQ. 1 )  DY1 = DY*0.25                                   
!      WY1 = DY1*2.0/(Y1-Y0)*WY                                          
!      Y01 =-DY1                                                         
!      Y11 = DY1                                                         
!                                                                       
!      CALL  TITLEP( TITLE,A )                                           
!      CALL  PLOT( XX,YY,-3 )                                            
!                                                                       
do 5 i=1,n
comp1(i) = 0.0d0
comp2(i) = 0.0d0
comp3(i) = 0.0d0
comp4(i) = 0.0d0
5 continue
!
!  ...  PLOT ORIGINAL DATA AND TREND  ...                               
!                                                                       
!      CALL  XYAXIS( 0.0,0.0,WX,WY,X0,X1,Y0,Y1,DX,DY,2,MESH )            
!      CALL  PLOTD( WORK,N,Y0,Y1,WX,WY,1,1,2 )                           
do 10 i=1,n
!xx   10 COMP1(I) = E(1,1,I)                                               
comp1(i) = e(1,1,i)
10 continue
!      CALL  PLOTD( Z,N,Y0,Y1,WX,WY,1,1,2 )                              
!      CALL  SYMBOL( WX-WC*18,WY+0.1,WC,'ORIGINAL AND TREND',0.0,18 )    
!                                                                       
!  ...  PLOT SEASONAL COMPONENT  ...                                    
!                                                                       
if( sorder .ne. 0) then
do 20 i=1,n
!xx   20 COMP2(I) = E(1,M1+M2+1,I)
comp2(i) = e(1,m1+m2+1,i)
20 continue
end if
!      CALL  MAXMIN( Z,N,ZMIN,ZMAX,DZ )                                  
!      DY2 = DY                                                          
!      IF( NNY .EQ. 1 )  DY2 = DY*0.5                                    
!      Y02 =-DY2                                                         
!      Y12 = DY2                                                         
!      NNY = 2                                                           
!  170 IF( Y02 .LE. ZMIN )  GO TO 180                                    
!      Y02 = Y02 - DY2                                                   
!      NNY = NNY + 1                                                     
!      GO TO 170                                                         
!  180 IF( Y12 .GE. ZMAX )  GO TO 190                                    
!      Y12 = Y12 + DY2                                                   
!      NNY = NNY + 1                                                     
!      GO TO 180                                                         
!  190 CONTINUE                                                          
!      WY2 = DY2*NNY*WY/(Y1-Y0)                                          
!      CALL  PLOT( 0.0,-WY2-1.5,-3 )                                     
!      CALL  XYAXIS( 0.0,0.0,WX,WY2,X0,X1,Y02,Y12,DX,DY2,2,MESH )        
!      CALL  PLOTD( Z,N,Y02,Y12,WX,WY2,1,1,2 )                           
!      CALL  SYMBOL( WX-WC*8,WY2+0.1,WC,'SEASONAL',0.0,8 )               
!                                                                       
!  ...  PLOT OBSERVATIONAL NOISE  ...                                   
!                                                                       

if( m4 .eq. 6 )  then
do 35 i=1,n
sum = 0.0d0
do 30 j=1,6
!xx   30 SUM = SUM + E(1,M123+J,N)*(TRADE(I,J)-TRADE(I,7))                 
sum = sum + e(1,m123+j,n)*(trade(i,j)-trade(i,7))
30 continue
!x   35 E(2,1,I) = SUM                                                    
e(2,1,i) = sum
35 continue
end if
if(m4 .eq. 1) then
do 37 i=1,n
tmp=trade(i,2)+trade(i,3)+trade(i,4)+trade(i,5)+trade(i,6)
!xx 37   e(2,1,i) = (trade(I,1)+trade(I,7)-0.4*tmp)*e(1,M123+1,n)
e(2,1,i) = (trade(i,1)+trade(i,7)-0.4*tmp)*e(1,m123+1,n)
37 continue
end if


!xx   40 IF( M5 .EQ. 0 )  GO TO 60                                         
if( m5 .eq. 0 )  go to 60
do 55 i=1,n
sum = 0.0d0
do 50 j=1,m5
!xx   50 SUM = SUM + E(1,M1234+J,N)*REG(I,J)                               
sum = sum + e(1,m1234+j,n)*reg(i,j)
50 continue
!xx   55 E(2,2,I) = SUM                                                    
e(2,2,i) = sum
55 continue
60 continue
!      CALL  PLOT( 0.0,-WY1-1.5,-3 )                                     
!      CALL  XYAXIS( 0.0,0.0,WX,WY1,X0,X1,Y01,Y11,DX,DY1,2,MESH )        
do 70 i=1,n
!c      Z(I) = WORK(I) - E(1,1,I)*ID(M1) - E(1,M1+1,I)*ID(M2)             
!c     *     - E(1,M12+1,I)*ID(M3) - E(2,1,I)*ID(M4) - E(2,2,I)*ID(M5)    
!c 70   COMP5(I)=Z(I)
!xx 70   COMP5(I) = Z(I) - E(1,1,I)*ID(M1) - E(1,M1+1,I)*ID(M2)
comp5(i) = z(i) - e(1,1,i)*id(m1) - e(1,m1+1,i)*id(m2)&
&- e(1,m12+1,i)*id(m3) - e(2,1,i)*id(m4) - e(2,2,i)*id(m5)
70 continue
!                                                                       
!      CALL  PLOTD( Z,N,Y01,Y11,WX,WY1,1,1,2 )                           
!      CALL  SYMBOL( WX-WC*5,WY1+0.1,WC,'NOISE',0.0,5 )                  
!                                                                       
if( m2 .eq. 0 )  go to 100
!     CALL  TITLEP( TITLE,A )                                           
!      CALL  PLOT( 16.0,WY1+WY2+3.0,-3 )                                 
!                                                                       
!  ...  PLOT AR PROCESS  ...                                            
!                                                                       
!      CALL  XYAXIS( 0.0,0.0,WX,WY2,X0,X1,Y02,Y12,DX,DY2,2,MESH )        
do 80 i=1,n
!xx   80 COMP3(I) = E(1,M1+1,I)                                            
comp3(i) = e(1,m1+1,i)
80 continue
!      CALL  PLOTD( Z,N,Y02,Y12,WX,WY2,1,1,2 )                           
!      CALL  SYMBOL( WX-WC*10,WY2+0.1,WC,'AR PROCESS',0.0,10 )           
!                                                                       
!  ...  PLOT DATA AND TREND PLUS AR  ...                                
!                                                                       
!      CALL  PLOT( 0.0,-8.0,-3 )                                         
!      CALL  XYAXIS( 0.0,0.0,WX,WY,X0,X1,Y0,Y1,DX,DY,2,MESH )            
!      CALL  PLOTD( WORK,N,Y0,Y1,WX,WY,1,1,2 )                           
!      DO 90 I=1,N                                                       
!   90 Z(I) = E(1,1,I) + E(1,M1+1,I)                                     
!      CALL  PLOTD( Z,N,Y0,Y1,WX,WY,1,1,2 )                              
!      CALL  SYMBOL(WX-WC*26,WY+0.1,WC,'ORIGINAL AND TREND PLUS AR',0.0, 
!     *             28 )                                                 
!                                                                       
100 if( m4 .eq. 0 )  go to 130
!      IF( M2 .NE. 0 )  CALL  PLOTI                                      
!      IF( M2 .NE. 0 )  CALL  TITLEP( TITLE,A )                          
!      IF( M2 .NE. 0 )  CALL  PLOT( XX,YY,-3 )                           
!      IF( M2 .EQ. 0 )  CALL  PLOT( 16.0,WY1+WY2+3.0,-3 )                
!                                                                       
!  ...  PLOT TRADING DAY EFFECT  ...                                    
!                                                                       
!      CALL  XYAXIS( 0.0,0.0,WX,WY2,X0,X1,Y02,Y12,DX,DY2,2,MESH )        
do 110 i=1,n
!xx  110 COMP4(I) = E(2,1,I)                                                   
comp4(i) = e(2,1,i)
110 continue
!      CALL  PLOTD( Z,N,Y02,Y12,WX,WY2,1,1,2 )                           
!      CALL  SYMBOL( WX-WC*18,WY2+0.1,WC,'TRADING DAY EFFECT',0.0,18 )   
!                                                                       
!  ...  PLOT SEASONAL PLUS TRADING DAY EFFECT  ...                      
!                                                                       
!      CALL  PLOT( 0.0,-WY2-2.0,-3 )                                     
!      CALL  XYAXIS( 0.0,0.0,WX,WY2,X0,X1,Y02,Y12,DX,DY2,2,MESH )        
!      DO 120 I=1,N                                                      
!  120 Z(I) = E(1,M12+1,I) + E(2,1,I)                                    
!      CALL  PLOTD( Z,N,Y02,Y12,WX,WY2,1,1,2 )                           
!      CALL  SYMBOL( WX-WC*32,WY2+0.1,WC,'SEASONAL PLUS TRADING DAY EFFEC
!     *T',0.0,32 )                                                       
130 continue
!      CALL  PLOTE                                                       
!      CALL  PLTCE                                                       
return
!xx  77  format(' ^*&$^$%#&^&*nvkdnfvafv')
end
!c      SUBROUTINE  REDATA( MT,DATA,ISW,IPR,X,N,TITLE,XM )                     
!xx      SUBROUTINE  REDATAD( DATA,ISW,IPR,X,N,XM )                     
subroutine  redatad( data,isw,x,n,xm )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE IS USED FOR THE LOADING OF ORIGINAL DATA,         
!     THE DATA IS LOADED THROUGH THE DEVICE SPECIFIED BY @MT@.          
!     EACH DATA SET IS COMPOSED OF TITLE, DATA LENGTH, DATA FORMAT      
!     AND ORIGINAL DATA.  IF ISW IS SET EQUAL TO 0, THE MEAN VALUE      
!     IS SUBTRACTED FROM THE ORIGINAL DATA.                             
!                                                                       
!       INPUTS:                                                         
!         MT:     INPUT DEVICE SPECIFICATION                            
!         ISW:    =0 MEAN VALUE IS SUBTRACTED FROM THE ORIGINAL DATA    
!                 =1 MEAN VALUE IS NOT SUBTRACTED                       
!         IPR:    =1 TO PRINT OUT DATA BY D-FORMAT                      
!                 =2 TO PRINT OUT DATA BY F-FORMAT                      
!         TITLE:  TITLE OF DATA                                         
!         N:      DATA LENGTH                                           
!         FORM:   INPUT DATA FORMAT SPECIFICATION                       
!         X(I) (I=1,N):  ORIGINAL DATA                                  
!                                                                       
!       OUTPUTS:                                                        
!         X:      ORIGINAL DATA                                         
!         N:      DATA LENGTH                                           
!         TITLE:  TITLE OF DATA                                         
!         XM:     MEAN                                                  
!                                                                       
!xx      IMPLICIT REAL*8( A-H,O-Z )                                        
!c      REAL * 4     FORM(20), TITLE(20)
!x      REAL*8       X(1) ,DATA(1)         
!xx      REAL*8       X(N) ,DATA(N)
integer isw, n
real(dp) data(n), x(n), xm
! local
integer i
real(dp) fn, s1, s2, s3, s4, xx
!                                                                       
!       LOADING OF TITLE, DATA LENGTH, FORMAT SPECIFICATION AND DATA    
!                                                                       
!      READ( MT,5 )      TITLE                                           
!      READ( MT,1 )      N                                               
!      READ( MT,5 )      FORM                                            
!      READ( MT,FORM )   (X(I),I=1,N)                                    
!                                                                       
!       ORIGINAL DATA PRINT OUT                                         
!                                                                       
do 53 i=1,n
!xx 53      X(I)=DATA(I)
x(i)=data(i)
53 continue

!      WRITE( 6,9 )     N , (FORM(I),I=1,17)                             
!      WRITE( 6,8 )     TITLE                                            
!      IF(IPR .EQ. 1)  WRITE( 6,3 )  (X(I),I=1,N)                        
!      IF(IPR .EQ. 2)  WRITE( 6,4 )  (X(I),I=1,N)                        
!                                                                       
!       SAMPLE MOMENTS COMPUTATION                                      
!   by Sato
!c	return



fn = n
s1 = 0.0d0
do 10     i=1,n
!xx   10 S1 = S1 + X(I)                                                    
s1 = s1 + x(i)
10 continue
xm = s1 / fn
s2 = 0.0d0
s3 = 0.0d0
s4 = 0.0d0
do 20  i=1,n
xx = x(i) - xm
s2 = s2 + xx**2
s3 = s3 + xx**3
!xx   20 S4 = S4 + XX**4                                                   
s4 = s4 + xx**4
20 continue
!                                                                       
s2 = s2 / fn
s3 = s3 / (fn*s2*dsqrt(s2))
s4 = s4 / (fn*s2*s2)
!                                                                       
!       PRINT OUT SAMPLE MOMENTS                                        
!                                                                       
!      WRITE( 6,7 )   XM, S2, S3, S4                                     
if( isw .eq. 1 )  return
!                                                                       
!       MEAN DELETION                                                   
!                                                                       
do 30  i=1,n
!xx   30 X(I) = X(I) - XM                                                  
x(i) = x(i) - xm
30 continue
!                                                                       
return
!                                                                       
!xx    1 FORMAT( 16I5 )                                                    
!xx    3 FORMAT( 1H ,10D13.5 )                                             
!xx    4 FORMAT( 1H ,10F10.4 )                                             

!xx    5 FORMAT ( 20A4 )                                                   
!xx    7 FORMAT( 1H0,'MEAN      =',D17.8,/,' VARIANCE  =',D17.8,/,         
!xx     *        ' SKEWNESS  =',D17.8,/,' KURTOSIS  =',D17.8 )             
!xx    8 FORMAT( 1H ,20A4 )                                                
!xx    9 FORMAT( 1H0,'<<  ORIGINAL DATA  X(I) (I=1,N)  >>',5X,'N =',I5,4X, 
!xx     *       'FORMAT =',17A4 )                                          
!                                                                       
end
subroutine  setfgh
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  SET F, G AND H MATRICES OF STATE SPACE MODEL  ...               
!                                                                       
!xx      IMPLICIT  REAL*8(A-H,O-Z)                                         
!xx      INTEGER   PERIOD, SORDER
! local
integer i
! common /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
! common /COMSM3/
real(dp) f1, f2, f3, a1, a2, a3, di, ui, tdf
!c      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,        
!c     *              NS, NI, MISING, IOUT, LL, N, NYEAR, nmonth,NPREDS,  
!c     *                    NPREDE, NPRED
!c      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
!c     *                   ,DI, TAU(3)                                   
common  /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
common  /comsm3/  f1(10), f2(10), f3(300), a1(10), a2(10),&
&a3(300), di, ui(3), tdf(7)
!                                                                       
di = 1.0d0
!                                                                       
if( m1 .eq. 0 )  go to 40
!xx      GO TO (10,20,30), M1                                              
if( m1 .eq. 1 ) go to 10
if( m1 .eq. 2 ) go to 20
if( m1 .eq. 3 ) go to 30
!                                                                       
10 f1(1) =  1.0d0
a1(1) =  1.0d0
go to 40
20 f1(1) = -1.0d0
f1(2) =  2.0d0
a1(1) =  2.0d0
a1(2) = -1.0d0
go to 40
30 f1(1) =  1.0d0
f1(2) = -3.0d0
f1(3) =  3.0d0
a1(1) =  3.0d0
a1(2) = -3.0d0
a1(3) =  1.0d0
!                                                                       
!c   40 IF( KSEA .EQ. 0 )  GO TO 160                                      
!c      IF( KSEA .EQ. 2 )  GO TO 130                                      
!c      IF( KSEA .EQ. -1 ) GO TO 80                                       
40 if( sorder .eq. 0 )  go to 160
if( sorder .eq. 2 )  go to 130
if( sorder .eq. -1 ) go to 80
do 70  i=1,m3
a3(i) = -1.0d0
!xx   70 F3(I) = -1.0D0
f3(i) = -1.0d0
70 continue
go to 150
!xx   80 DO 90 I=1,M3                                                      
!xx      A3(I) = 0.0D0                                                     
!xx   90 F3(I) = 0.0D0
80 continue
a3(1:m3) = 0.0d0
f3(1:m3) = 0.0d0
a3(m3)= 1.0d0
f3(1) = 1.0d0
go to 150
!c  130 DO 140  I=1,ISEA-1                                                
130 do 140  i=1,period-1
f3(i) = -i
a3(i) = -i - 1
!c      F3(ISEA+I-1) = I - ISEA - 1                                       
!c  140 A3(ISEA+I-1) = I - ISEA                                           
f3(period+i-1) = i - period - 1
!xx  140 A3(PERIOD+I-1) = I - PERIOD
a3(period+i-1) = i - period
140 continue
150 continue
160 continue
!                                                                       
return
!                                                                       
!                                                                       
end
subroutine  smoth3( z,r,s,t,trade,imis,n,lm1,f,ilkf,ismt )
  use timsac_kinds, only: dp
  implicit none
!c      SUBROUTINE  SMOTH3( Z,R,S,T,TRADE,REG,EPRED,VPRED,IMIS,N,LM1,
!c     *                    F,ILKF,ISMT )
!C     *                    MJ,F,ILKF,ISMT )                              
!                                                                       
!     ...  INFORMATION SQUARE ROOT FILTER & SMOOTHER  ...               
!          FOR SEASONAL ADJUSTMENT                                      
!                                                                       
!xx      IMPLICIT  REAL*8  ( A-H,O-Z )                                     
!c      REAL*8     Z(N), T(L,LM1,N), REG(N,M5), EPRED(N), VPRED(N), TRADE
!c      DIMENSION  S(MJ,MJ), R(MJ,MJ), IMIS(N), TRADE(N,7)                
!c      DIMENSION  D(100), WI(10,10), X(40), TT(40), E(40)
!xx      REAL*8     Z(N), T(L,LM1,N), REG(N,M5)
!xx      DIMENSION  S(M+1,M+1), R(LM1+1,LM1), IMIS(N), TRADE(N,7)            
!xx      DIMENSION  D(LM1+1), WI(3,3), X(LM1), TT(L), E(M)                
!xx      INTEGER    PERIOD, SORDER
integer n, lm1, ilkf, ismt, imis(n)
! common /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
real(dp) z(n), r(lm1+1,lm1), s(m+1,m+1), t(l,lm1,n),&
&trade(n,7), f
! local
integer i, ii, iii, id, j, jj, l1, l12, l123, l1234, ll2, lm, m12,&
&m123, m45, nnn
real(dp) reg(n,m5), d(lm1+1), wi(3,3), x(lm1), tt(l),&
&e(m), pai, sdet, sum, tmp, tmpp
! common /COMSM3/
real(dp) f1, f2, f3, a1, a2, a3, di, ui, tdf
! common /CMFUNC/
real(dp) djacob, fc, sig2, aic, fi, sig2i, aici, gi, gc
!c      COMMON    /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,        
!c     *             NS, NI, MISING, IOUT, LL, NN, NYEAR,nmonth, NPREDS, 
!c     *                    NPREDE, IPRED
!cx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
!      COMMON    /CMFUNC/  DJACOB, FF, SIG2, AIC                        
!c      COMMON    /CMFUNC/  DJACOB,FF,
!c     *            SIG2,AIC,FI,SIG2I,AICI,GI(20),G(20) 
!cx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
common  /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
common  /comsm3/  f1(10), f2(10), f3(300), a1(10), a2(10),&
&a3(300), di, ui(3), tdf(7)
common  /cmfunc/  djacob,fc,sig2,aic,fi,sig2i,aici,gi(200),gc(200)

data  pai/3.1415926535d0/
!                                                                       
!c      M    = M1 + M2 + M3 + M4 + M5                                     
m12  = m1 + m2
m123 = m1 + m2 + m3
m45 = m4 + m5
lm   = l + m
l1   = l + m1
l12  = l + m1 + m2
l123 = l + m1 + m2 + m3
l1234= l + m1 + m2 + m3 + m4
ll2  = id(m1) + id(m2)
sig2 = 0.0d0
sdet = 0.0d0
nnn = 0
!
do 500  ii=1,n
!                                                                       
!     ...  TIME  UPDATE (PREDICTION)  ...                               
!                                                                       
if( m1 .eq. 0 )  go to 140
if(m1.eq.1)  go to 115
!xx      DO 110  J=2,M1                                                    
do 111  j=2,m1
do 110  i=1,j
!xx  110 R(L+I,L+J-1) = S(I,J)                                             
r(l+i,l+j-1) = s(i,j)
110 continue
111 continue
115 do 130  i=1,m1
sum = 0.0d0
do 120  j=i,m1
!xx  120 SUM = SUM + S(I,J)*F1(J)                                          
sum = sum + s(i,j)*f1(j)
120 continue
!xx  130 R(L+I,L1) = SUM                                                   
r(l+i,l1) = sum
130 continue
!                                                                       
140 if( m2 .eq. 0 )  go to 5
if(m2 .eq. 1)  go to 155
!xx      DO 150  J=2,M2                                                    
do 151  j=2,m2
jj = m1 + j
do 150  i=1,jj
!xx  150 R(L+I,L1+J-1) = S(I,M1+J)                                         
r(l+i,l1+j-1) = s(i,m1+j)
150 continue
151 continue
155 do 170  i=1,m2
sum = 0.0d0
do 160  j=i,m2
!xx  160 SUM = SUM + S(M1+I,M1+J)*F2(J)                                    
sum = sum + s(m1+i,m1+j)*f2(j)
160 continue
!xx  170 R(L1+I,L12) = SUM                                                 
r(l1+i,l12) = sum
170 continue
do 190  i=1,m1
sum = 0.0d0
do 180  j=1,m2
!xx  180 SUM = SUM + S(I,M1+J)*F2(J)                                       
sum = sum + s(i,m1+j)*f2(j)
180 continue
!xx  190 R(L+I,L12) = SUM                                                  
r(l+i,l12) = sum
190 continue
!xx  195 CONTINUE                                                          
!                                                                       
5 if( m3 .eq. 0 )  go to 55
!xx      DO 10  J=2,M3                                                     
do 11  j=2,m3
jj = m12 + j
do 10  i=1,jj
!xx   10 R(L+I,L12+J-1) = S(I,M12+J)                                       
r(l+i,l12+j-1) = s(i,m12+j)
10 continue
11 continue
do 30  i=1,m3
sum = 0.0d0
do 20  j=i,m3
!xx   20 SUM = SUM + S(M12+I,M12+J)*F3(J)                                  
sum = sum + s(m12+i,m12+j)*f3(j)
20 continue
!xx   30 R(L12+I,L123) = SUM                                               
r(l12+i,l123) = sum
30 continue
!
do 50  i=1,m12
sum = 0.0d0
do 40  j=1,m3
!xx   40 SUM = SUM + S(I,M12+J)*F3(J)                                      
sum = sum + s(i,m12+j)*f3(j)
40 continue
!xx   50 R(L+I,L123) = SUM                                                 
r(l+i,l123) = sum
50 continue
!                                                                       
55 if( m45 .eq. 0 )  go to 70
!xx      DO 60  J=1,M45
do 61  j=1,m45
jj = m123 + j
do 60  i=1,jj
!xx   60 R(L+I,L123+J) = S(I,M123+J)
r(l+i,l123+j) = s(i,m123+j)
60 continue
61 continue
70 continue
!                                                                       
if( m1 .eq. 0 )  go to 215
do 210  i=1,m1
!xx  210 R(L+I,1) = -R(L+I,L+1)                                            
r(l+i,1) = -r(l+i,l+1)
210 continue
215 if( m2 .eq. 0 ) go to 225
do 220  i=1,m12
!xx  220 R(L+I,LL2) = -R(L+I,L+M1+1)
r(l+i,ll2) = -r(l+i,l+m1+1)
220 continue
!                                                                       
225 if( m3 .eq. 0 )  go to 235
do 230  i=1,m123
!xx  230 R(L+I,L) = -R(L+I,L+M12+1)                                        
r(l+i,l) = -r(l+i,l+m12+1)
230 continue
!                                                                       
!xx  235 DO 250  I=1,L                                                     
!xx      DO 250  J=1,LM+1                                                  
!xx  250 R(I,J) = 0.0D0
235 continue
r(1:l,1:lm+1) = 0.0d0
do 240  i=1,m
!xx  240 R(L+I,L+M+1) = S(I,M+1)                                           
r(l+i,l+m+1) = s(i,m+1)
240 continue
!                                                                       
do 260 i=1,l
!xx  260 R(I,I) = UI(I)                                                    
r(i,i) = ui(i)
260 continue
!                                                                       
!c      CALL  HUSHL7( R,D,MJ,L+M+1,L12+1,L+M123 )                         
call  hushl7( r,d,lm1+1,lm1,l12+1,l+m123 )
!                                                                       
!xx      DO 300  I=1,L                                                     
do 301  i=1,l
do 300  j=1,lm+1
!xx  300 T(I,J,II) = R(I,J)
t(i,j,ii) = r(i,j)
300 continue
301 continue
!                                                                       
!     ...  MEASUREMENT  UPDATE (FILTERING)  ...                         
!                                                                       
!    modified at 96.10 by S.S.
! 400  continue
!xx 400    IF( IMIS(II) .EQ. 1 )  GO TO 420
if( imis(ii) .eq. 1 )  go to 420
!
nnn=nnn+1
!     WRITE(6,997) M,UI(1)                                              
!xx      DO 410  J=1,M+1                                                   
do 411  j=1,m+1
s(m+1,j) = 0.0d0
do 410  i=1,j
!xx  410 S(I,J) = R(L+I,L+J)                                               
s(i,j) = r(l+i,l+j)
410 continue
411 continue
!                                                                       
s(m+1,1)    = di
s(m+1,m1+1) = di
s(m+1,m12+1)= di
s(m+1,m+1)  = di*z(ii)
if( m4 .eq. 6 )  then
do 414  i=1,6
jj = ii
if( ilkf .eq. 0 )  jj = n - ii + 1
!xx  414 S(M+1,M123+I) = DI*(TRADE(JJ,I) - TRADE(JJ,7))                    
s(m+1,m123+i) = di*(trade(jj,i) - trade(jj,7))
414 continue
end if

if(m4 .eq. 1) then

jj = ii
if( ilkf .eq. 0 )  jj = n - ii + 1
tmpp=trade(jj,2)+trade(jj,3)
tmpp=tmpp+trade(jj,4)+trade(jj,5)+trade(jj,6)
s(m+1,m123+1)=di*(trade(jj,1)+trade(jj,7)-0.4*tmpp)
end if

!xx  415 IF( M5 .EQ. 0 )  GO TO 418                                        
if( m5 .eq. 0 )  go to 418
do 417  i=1,m5
!xx  417 S(M+1,M123+M4+I) = REG(I,II)                                      
s(m+1,m123+m4+i) = reg(i,ii)
417 continue
418 continue
!                                                                       
!     CALL  PRINT2( S,M+1,M+1,MJ,1 )                                    
!c      CALL  HUSHL4( S,MJ,M+1,M+1,1,0 )                                  
call  hushl4( s,m+1,m+1,m+1,1,0 )
!     CALL  PRINT2( S,M+1,M+1,MJ,1 )                                    
!                                                                       
!     ...  PREDICTION ERROR & DETERMINANT  ...                          
!                                                                       
!c      IF( IPRED.EQ.1 .AND. II.GE.NPREDS )  GO TO 430                    
if(ilkf .eq. 0)  go to 500
sig2 = sig2 + s(m+1,m+1)**2
do 440  i=1,m
!xx  440 SDET = SDET + DLOG( S(I,I)**2 ) - DLOG( R(L+I,L+I)**2 )           
sdet = sdet + dlog( s(i,i)**2 ) - dlog( r(l+i,l+i)**2 )
440 continue
go to 500
!                                                                       
!  ...  LONG RANGE PREDICTION  ...                                      
!                                                                       
!c  430 CALL  RECOEF( R,LM,LM,MJ,X )                                      
!xx  430 CALL  RECOEF( R,LM,LM,LM1+1,X )
call  recoef( r,lm,lm,lm1+1,x )
sum = id(m1)*x(l+1) + id(m2)*x(l1+1) + id(m3)*x(l12+1)
if( m4 .eq. 6 )  then
do 431 j=1,6
!xx  431 SUM = SUM + (TRADE(II,J)-TRADE(II,7))*X(L123+J)
sum = sum + (trade(ii,j)-trade(ii,7))*x(l123+j)
431 continue
end if
if(m4 .eq. 1) then
tmp=trade(ii,2)+trade(ii,3)+trade(ii,4)+trade(ii,5)+trade(ii,6)
sum=sum+(trade(ii,1)+trade(ii,7)-0.4*tmp)*x(l123+1)
end if
!c  432 IF( M5 .EQ. 0 )  GO TO 434                                        
!c      DO 433 J=1,M5                                                     
!c  433 SUM = SUM + REG(II,J)*X(L1234+J)                                  
!c  434 EPRED(II) = SUM                                                   
!c      SUM = 1.0D0                                                       
!c      DO 445 I=1,M                                                      
!c  445 SUM = SUM * S(I,I)/R(L+I,L+I)                                     
!c      VPRED(II) = SUM**2                                                
!cc      WRITE(6,977)  II, EPRED(II), VPRED(II)                            
!c  977 FORMAT( 1H ,I5,2F15.7 )                                           
!                                                                       
!  ...  MISSING OBSERVATION  ...                                        
!                                                                       
!xx  420 DO 425  I=1,M
420 do 426  i=1,m
do 425  j=1,m+1
!xx  425 S(I,J) = R(L+I,L+J)
s(i,j) = r(l+i,l+j)
425 continue
426 continue
500 continue
!                                                                       
!     ... LIKELIHOOD COMPUTATION  ...                                   
!                                                                       

!   99/8/12
if(ilkf .eq. 0)  go to 600

sig2 = sig2/(2.0d0*nnn)
f = -0.5d0*(nnn*( dlog(pai*2.d0) + dlog(sig2) + 1.0d0 ) + sdet)
!c      FF = F                                                            
fc = f
!    modified 97/10/17
aic = -2.0d0*(f+djacob) + 2.0d0*(m2 + l + 1 + m)
!xx  901 FORMAT(1H ,4D16.9)                                                
! 900 FORMAT(1H ,'F =',D15.8,5X,'SIG2 =',D15.8,5X,'DET =',D15.8,        
!    1     5X,'FF =',D15.8,5X,'SD =',D15.8,5X,'LOG V =',D15.8 )         
!xx  900 FORMAT(1H ,3D20.10)                                               
!                                                                       
!                                                                       
!                                                                       
!     ...  SMOOTHING  ...                                               
!                                                                       
600 if(ismt .eq.  0) return
!c      CALL  RECOEF( S,M,M,MJ,X )
call  recoef( s,m,m,m+1,x )
!                                                                       
do 610 i=1,m
!xx  610 E(I) = X(I)
e(i) = x(i)
610 continue
!                                                                       
do 1200  iii=1,n-1
!                                                                       
ii = n-iii+1
!                                                                       
!  ...  WI = INVERSE OF T(I,J)  ...                                     
!                                                                       
wi(1,1) = 1.0d0/t(1,1,ii)
if( l .eq. 1 )  go to 615
wi(2,2) = 1.0d0/t(2,2,ii)
wi(1,2) = -t(1,2,ii)*wi(2,2)*wi(1,1)
if( l .eq. 2 )  go to 615
wi(3,3) = 1.0d0/t(3,3,ii)
wi(2,3) = -t(2,3,ii)*wi(2,2)*wi(3,3)
wi(1,3) = -(t(1,2,ii)*wi(2,3) + t(1,3,ii)*wi(3,3))*wi(1,1)
615 do 630  i=1,l
sum = t(i,lm1,ii)
do 620  j=1,m
!xx  620 SUM = SUM - T(I,L+J,II)*X(J)
sum = sum - t(i,l+j,ii)*x(j)
620 continue
!xx  630 TT(I) = SUM
tt(i) = sum
630 continue
!                                                                       
do 650  i=1,l
sum = 0.0d0
do 640  j=i,l
!xx  640 SUM = SUM + WI(I,J)*TT(J)                                         
sum = sum + wi(i,j)*tt(j)
640 continue
!xx  650 D(I) = SUM
d(i) = sum
650 continue
x(1) = x(1) - d(1)
if(m1.gt.0 .and. m2.gt.0)  x(m1+1) = x(m1+1) - d(2)
x(m12+1) = x(m12+1) - d(l)
do 660  i=1,m
!xx  660 D(I) = X(I)
d(i) = x(i)
660 continue
!                                                                       
if( m1 .eq. 0 )  go to 715
do 710  i=1,m1
!xx  710 X(I) = F1(I)*D(M1)                                                
x(i) = f1(i)*d(m1)
710 continue
715 if(m1.le.1)  go to  730
do 720  i=2,m1
!xx  720 X(I) = X(I) + D(I-1)                                              
x(i) = x(i) + d(i-1)
720 continue
!                                                                       
730 if( m2 .eq. 0 )  go to 760
do 740  i=1,m2
!xx  740 X(M1+I) = F2(I)*D(M12)
x(m1+i) = f2(i)*d(m12)
740 continue
if(m2 .le. 1)  go to 760
do 750  i=2,m2
!xx  750 X(M1+I) = X(M1+I) + D(M1+I-1)                                     
x(m1+i) = x(m1+i) + d(m1+i-1)
750 continue
!                                                                       
760 if( m3 .le. 0 )  go to 775
do 770  i=1,m3
!xx  770 X(M12+I) = F3(I)*D(M123)                                          
x(m12+i) = f3(i)*d(m123)
770 continue
775 if( m3 .le. 1 )  go to 785
do 780  i=2,m3
!xx  780 X(M12+I) = X(M12+I) + D(M12+I-1)
x(m12+i) = x(m12+i) + d(m12+i-1)
780 continue
!                                                                       
785 continue
!xx 2000 CONTINUE                                                          
do 790 i=1,m
t(1,i,ii) = e(i)
!xx  790 E(I) = X(I)                                                       
e(i) = x(i)
790 continue
!                                                                       
1200 continue
!                                                                       
do 1210 i=1,m
!xx 1210 T(1,I,1) = E(I)                                                   
t(1,i,1) = e(i)
1210 continue
if( m4 .eq. 6 ) then
sum = 0.0d0
do 1220 i=1,6
tdf(i) = x(m123+i)
!xx 1220 SUM = SUM + TDF(I)                                                
sum = sum + tdf(i)
1220 continue
tdf(7) = -sum
end if
if( m4 .eq. 1 ) then
tdf(1) = x(m123+1)
tdf(7) = x(m123+1)
do 1222 i=2,6
!xx 1222    TDF(i) = -0.4*TDF(1)
tdf(i) = -0.4*tdf(1)
1222 continue
end if
!                                                                       
return
end
!xx      SUBROUTINE  SPARAM0( N,IPAR,NIP,para,NPA )
subroutine  sparam0( ipar,nip )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  Set or read control parameters ...                              
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION  IPAR(NIP), para(NPA)
!xx      INTEGER    PERIOD, SORDER, TRADE
integer nip, ipar(nip)
! local
integer id, trade
! common /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
! common /CCC/
integer isw, ismt, idif, log, mesh
common    /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
!xx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
common    /ccc/     isw, ismt, idif, log, mesh
!                                                                       
!  ...  SET DEFAULT VALUES  ...                                         
!                                                                       
m1 = ipar(1)
m2 = ipar(2)
m4 = 0
m5 = 0
log  = ipar(5)
!xx      IPR  = 7                                            
isw  = 1
idif = ipar(7)
if(idif .gt. 2) idif=1
if(idif .lt. 1) idif=1
mesh = 1
period = ipar(3)
sorder = ipar(4)
trade  = ipar(6)
if(trade .eq. 1) trade = 7
if(trade .ge. 1)  nyear=ipar(8)
if(trade .ge. 1)  nmonth=ipar(9)
m3 = (period - 1)*sorder
if(sorder.eq.-1)  m3 = period
if(trade .ge. 1)  m4 = trade-1
m = m1 + m2 + m3 + m4 + m5
l = id(m1) + id(m2) + id(m3)
l = max(2,l)
!                                                                       
return
end
!c      SUBROUTINE  SPARAM( MT,A,IPAR,para,iopt )                                        
!xx      SUBROUTINE  SPARAM( N,A,NA,IPAR,NIP,para,NPA,iopt )
subroutine  sparam( a,na,para,npa,iopt )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  Set or read control parameters ...                              
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      REAL*8     TAU2(3), A(40), PAC(10), F1,F2,F3,A1,AR,A3,para(26)
!c      REAL*8     ARCC(10)
!c      INTEGER   PERIOD, SORDER, BSPAN, OUTLIR, TRADE, YEAR, PRED        
!c      INTEGER   PREDS, PREDE ,IPAR(11),month
!c      DIMENSION  TAU2(3), PAC(10), ARCC(10)
!xx      DIMENSION  A(NA), IPAR(NIP), para(NPA)
!xx      DIMENSION  TAU2(3), PAC(M2), ARCC(M2)
!xx      INTEGER    PERIOD, SORDER
integer na, npa, iopt
real(dp) a(na), para(npa)
! common /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
! local
integer i
real(dp) tau2(3), pac(m2), arcc(m2)
! common /COMSM3/
real(dp) f1, f2, f3, a1, a2, a3, di, ui, tdf
! common /CCC/
integer isw, ismt, idif, log, mesh
!c     *              BSPAN, ISPAN, MISING, OUTLIR, LL, N, YEAR, month,
!c     *                    PREDS, PREDE, PRED
!c      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), AR(10), A3(30)    
!cx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)    
!xx      COMMON    /CCC/     ISW, IPR, ISMT, IDIF, LOG, MESH              
common    /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
common    /comsm3/  f1(10), f2(10), f3(300), a1(10), a2(10),&
&a3(300), di, ui(3), tdf(7)
common    /ccc/     isw, ismt, idif, log, mesh
!      NAMELIST  /PARAM/  M1, M2, M5, PERIOD, SORDER, TRADE, MT, BSPAN,  
!     *                   ISPAN, MISING, OUTLIR, TAU2, PAC, IPR, IDIF,   
!     *                   LOG, YEAR, PRED, MESH                          
!                                                                       
!  ...  SET DEFAULT VALUES  ...                                         
!                                                                       
!c      M1 = IPAR(1)       
!c      M2 = IPAR(2)
!c      m4 = 0
!c      M5 = 0                        
!c      MT = 1                                                            
!c      LOG  = IPAR(5)                            
!c      IPR  = 7                                            
!c      ISW  = 1                                                          
!c      IDIF = IPAR(7)
!c      if(idif .gt. 2) idif=1
!c      if(idif .lt. 1) idif=1
!c      MESH = 1                                                          
!c      PERIOD = IPAR(3)                               
!c      SORDER = IPAR(4)                                                 
!c	TRADE  = IPAR(6)                
!c      if(TRADE .EQ. 1) TRADE = 7
!c      IF(TRADE .GE. 1)  YEAR=IPAR(10)
!c      IF(TRADE .GE. 1)  month=IPAR(11)
!      BSPAN  = IPAR(8)                                                 
!      ISPAN  = IPAR(9)          
!      if(n .lt. ispan) ispan = n
!c      bspan = n
!c      ispan = n
!c      if(iopt .lt. 0) ispan = ipar(9)
!c      PRED = 0                                                          
!c      MISING = 0                                                        
!c      OUTLIR = 0                                                        
tau2(1) = 0.005d0
tau2(2) = 0.800d0

! 99/8/12    
tau2(3) = 1.0d-3
!      TAU2(3) = 1.0D-5    
if(m2 .eq. 0) tau2(2) = 0.001d0

!xx      do 18 i=1,7
!xx 18      TDF(i) = 0.0D00
tdf(1:7) = 0.0d00
!c      DO 20  I=1,10                                                     
do 20  i=1,m2
!xx   20 PAC(I) = 0.88D0*(-0.6D0)**(I-1)                                   
pac(i) = 0.88d0*(-0.6d0)**(i-1)
20 continue

if(iopt .lt. 0) then
do 21 i=1,3
tau2(i)=para(i) - 0.0001
if(tau2(i) .ge. 1.0d0) tau2(i)=1.0d0 - 0.1d-20
if(tau2(i) .le. 0.0d0) tau2(i)=0.0d0 + 0.1d-20
21 continue

!            write(*,*) tau2(1),tau2(2),tau2(3)
if(m2 .gt. 0) then
do 22 i=1,m2
!xx 22            arcc(i) = para(3+i)
arcc(i) = para(3+i)
22 continue
call parcor( arcc,m2,pac )
!            write(*,*) PAC(1),PAC(2)
end if
end if
!                                                                       
!  ...  READ IN CONTROL PARAMETERS  ...                                 
!                                                                       
!      READ(5,PARAM)                                                     
!                                                                       
!c      M3 = (PERIOD - 1)*SORDER                                          
!c      IF(SORDER.EQ.-1)  M3 = PERIOD                                     
!c      IF(TRADE .GE. 1)  M4 = TRADE-1                                          
!c      M = M1 + M2 + M3 + M4 + M5                                        
!c      L = ID(M1) + ID(M2) + ID(M3)                                      
!c      LL= ID(M1) + ID(M2) + ID(M3) + ID(M4) + ID(M5)                    
do 30 i=1,l
!xx   30 A(I) = DASIN( TAU2(I)*2.0D0 - 1.0D0 )                             
a(i) = dasin( tau2(i)*2.0d0 - 1.0d0 )
30 continue
!      WRITE(6,610)                                                      
!      WRITE(6,600) M1,M2,M3,M4,M5,M,L,PERIOD,SORDER,BSPAN,ISPAN,MISING, 
!     *           OUTLIR, LL                                             
if(m2.eq.0)  return
!                                                                       
do 40 i=1,m2
!xx   40 A(L+I) = DASIN( PAC(I)/0.90D0 )
a(l+i) = dasin( pac(i)/0.90d0 )
40 continue
!                                                                       
return
!xx  600 FORMAT( 10X,'M1     =',I3,/,10X,'M2     =',I3,/,10X,'M3     =',   
!xx     *   I3,/,10X,'M4     =',I3,/,10X,'M5     =',I3,/,10X,'M      =',   
!xx     *   I3,/,10X,'L      =',I3,/,10X,'PERIOD =',I3,/,10X,'SORDER =',   
!xx     *   I3,/,10X,'BSPAN  =',I3,/,10X,'ISPAN  =',I3,/,10X,'MISING =',   
!xx     *   I3,/,10X,'OUTLIR =',I3,/,10X,'LL     =',I3 )                   
!xx  610 FORMAT( 5X,'---  PROGRAM  DECOMP  ---' )                          
end
subroutine  state( x,a,k )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!  ...  TRANSFORMATION OF STATE VECTOR FOR TIME REVERSED MODEL  ...     
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!c      DIMENSION  X(K), A(K), Y(30)                                      
!xx      DIMENSION  X(K), A(K), Y(K)
integer k
real(dp) x(k), a(k)
! local
integer i, j
real(dp) y(k), sum
!                                              
y(1) = a(1)*x(1)
y(2:k) = 0.0d0
if( k .eq. 0 )  return
sum = 0.0d0
!xx      DO 20  I=1,K
do 21  i=1,k
sum = a(i)*x(1)
if( i .lt. k )  sum = sum + x(i+1)
if( i .eq. 1 )  go to 20
do 10  j=1,i-1
!xx   10 SUM = SUM + A(J)*Y(I-J)
sum = sum + a(j)*y(i-j)
10 continue
20 y(i) = sum
21 continue
!                                                                       
x(1) = y(1)
if( k .eq. 1 ) return
do 40  i=2,k
sum = 0.0d0
do 30  j=i,k
!xx   30 SUM = SUM + A(J)*Y(J-I+2)                                         
sum = sum + a(j)*y(j-i+2)
30 continue
!xx   40 X(I) = SUM
x(i) = sum
40 continue
!                                                                       
return
end

subroutine  trade( jsyear,nmonth,n,tday )
  use timsac_kinds, only: dp
  implicit none

!xx      IMPLICIT REAL*8(A-H,O-Z)                                          

!
!  ...  This subroutine computes the number of days of the week
!       in each month, Nov.1981
!    modified at '96 by S.S.
!    This subroutine should not be used after 2099.
!
!xx      DIMENSION  TDAY(N,7), IX(12)
integer jsyear, nmonth, n
real(dp) tday(n,7)
! local
integer i, ie, ii, i0, i1, i2, j, jj, js
integer ix(12)
data   ix  /3,0,3,2,3,2,3,3,2,3,2,3/
!
!      open(1,file='tmp.dat')
js = jsyear - 1900
!c      I0 = MOD( JS+(JS-1)/4,7 ) + 1
i2 = mod( js+(js-1)/4,7 ) + 1
jj = 2-nmonth
ii = 0
5 ii = ii + 1
i1 = ii + js - 1
ix(2) = 0
if( mod(i1,4) .eq. 0 )  ix(2) = 1
if( mod(i1+1900,100).eq.0 )  ix(2) = 0
if( mod(i1+1900,400).eq.0 )  ix(2) = 1
do 30  j=1,12
do 10  i=1,7
!xx 10    if( jj.gt.0 ) TDAY(JJ,I) = 4.0
if( jj.gt.0 ) tday(jj,i) = 4.0
10 continue
!
ie = ix(j)
if( ie .eq. 0 )  go to 28
i0 = i2
do 20  i=1,ie
i2 = i0 + i
if( i2 .gt. 7 ) i2 = i2 - 7
!xx 20    if( jj.gt.0 ) TDAY(JJ,I2) = 5.0
if( jj.gt.0 ) tday(jj,i2) = 5.0
20 continue
!c      I0 = I2

!       if(jj .gt.0) WRITE(1,*)  (TDAY(JJ,I),I=1,7)
28 jj = jj + 1
if( jj .gt.n ) then

!      do 29 i9=1,n
! 29      write(1,*) (tDAY(i9,j9),j9=1,7)

return
end if
30 continue
go to 5
!
end

subroutine  trade2( jsyear,nquart,n,tday )
  use timsac_kinds, only: dp
  implicit none

!xx      IMPLICIT REAL*8(A-H,O-Z)                                          

!
!  ...  This subroutine computes the number of days of the week
!       in each quarter.
!    modified at '96 by S.S.
!    This subroutine should not be used after 2099.
!
!xx      DIMENSION  TDAY(N,7), IX(4)
integer jsyear, nquart, n
real(dp) tday(n,7)
! local
integer i, ie, ii, i0, i1, i2, j, jj, js
integer ix(4)
data   ix  /6,7,8,8/
!
!      write(6,*) nquart, jsyear
js = jsyear - 1900
!c      I0 = MOD( JS+(JS-1)/4,7 ) + 1
i2 = mod( js+(js-1)/4,7 ) + 1
jj = 2-nquart
ii = 0
5 ii = ii + 1
i1 = ii + js - 1
ix(1) = 6
if( mod(i1,4) .eq. 0 )  ix(1) = 7
if( mod(i1+1900,100).eq.0 )  ix(1) = 6
if( mod(i1+1900,400).eq.0 )  ix(1) = 7
do 30  j=1,4
do 10  i=1,7
!xx 10    if( jj.gt.0 ) TDAY(JJ,I) = 12.0
if( jj.gt.0 ) tday(jj,i) = 12.0
10 continue
!
ie = ix(j)
!      write(6,*) IE, jj
if( ie .eq. 0 )  go to 28
i0 = i2
do 20  i=1,ie
i2 = i0 + i
if( i2 .gt. 7 ) i2 = i2- 7
if( i2 .gt. 7 ) i2 = i2- 7
!      write(6,*) i2, tday(jj,i2)
!xx 20    if( jj.gt.0 ) TDAY(JJ,I2) = TDAY(JJ,I2) + 1.0
if( jj.gt.0 ) tday(jj,i2) = tday(jj,i2) + 1.0
20 continue
!c      I0 = I2
!      WRITE(6,*)  (TDAY(JJ,I),I=1,7)
28 jj = jj + 1
if( jj .gt.n ) return
30 continue
go to 5
!
end

!c      subroutine trpar( a,para )
subroutine trpar( a,na,para,npa )
  use timsac_kinds, only: dp
  implicit none
!      SUBROUTINE  TITLEP( TITLE,A )                                     
!                                                                       
!  ...  PLOT DATA ID AND ESTIMATED PARAMETERS  ...                      
!                                                                       
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!      REAL*4   FM1, FM2, FM3, FM4, FM5, FAIC,FSIG2,TAU1,TAU2,TAU3       
!      REAL*4     TITLE(20), DAY(2), TIME(3)                             
!c      DIMENSION  A(40),para(26),a2(10),atmp(10)
!xx      DIMENSION  A(na),para(npa),ar(M2),atmp(M2)
!xx      INTEGER    PERIOD, SORDER
integer na, npa
real(dp) a(na), para(npa)
! common /COMSM2/
integer m1, m2, m3, m4, m5, m, l, period, sorder, nyear, nmonth
! local
integer i
real(dp) ar(m2), atmp(m2), tau1, tau2, tau3
! common /COMSM3/
real(dp) f1, f2, f3, a1, a2, a3, di, ui, tdf
! common /CMFUNC/
real(dp) djacob, fc, sig2, aic, fi, sig2i, aici, gi, gc
!c      COMMON  /COMSM2/  M1, M2, M3, M4, M5, M, L, ISEA, KSEA,                  
!c     *         NS,NI,MISING,IOUT,LL,N,NYEAR,nmonth,NPS,NPE,NPRED    
!      COMMON  /CMFUNC/  DJACOB, F, SIG2, AIC                           
!c      COMMON    /CMFUNC/  DJACOB,F,SIG2,AIC,FI,SIG2I,AICI,GI(20),G(20) 
!c      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), AR(10), A3(30)    
!cx      COMMON    /CMFUNC/  DJACOB,FC,SIG2,AIC,FI,SIG2I,AICI,GI(20),GC(20)
!cx      COMMON    /COMSM3/  F1(10), F2(10), F3(30), A1(10), A2(10), A3(30)
common  /comsm2/  m1, m2, m3, m4, m5, m, l, period, sorder,&
&nyear, nmonth
common  /cmfunc/  djacob,fc,sig2,aic,fi,sig2i,aici,gi(200),gc(200)
common  /comsm3/  f1(10), f2(10), f3(300), a1(10), a2(10),&
&a3(300), di, ui(3), tdf(7)
!                                                                       
tau1 = 0.0d00
tau2 = 0.0d00
tau3 = 0.0d00


!      FM1 = M1                                                          
!      FM2 = M2                                                          
!      FM3 = M3                                                          
!      FM4 = M4                                                          
!      FM5 = M5
para(1) = aic
!c      para(2) = F                                                        
para(2) = fc
para(3) = sig2
!      CALL  DATE( DAY )                                                 
!      CALL  CLOCK( TIME,1 )                                             
!     TAU1 = DEXP( A(1) )/(1.0D0 + DEXP(A(1)))                          
!     IF(L.GE.2)  TAU2 = DEXP( A(2) )/(1.0D0 + DEXP( A(2) ))            
!     IF(L.GE.3)  TAU3 = DEXP( A(3) )/(1.0D0 + DEXP( A(3) ))            
tau1 = 0.5d0*(1.0d0 + dsin( a(1) )) + 0.0001

if(l.ge.2)  tau2 = 0.5d0*(1.0d0 + dsin( a(2) )) + 0.0001
!     IF(L.EQ.3)  TAU2 = DEXP( A(2) )                                   
if(l.ge.3)  tau3 = 0.5d0*(1.0d0 + dsin( a(3) )) + 0.0001
!     TAU1 = A(1)**2                                                    
!     TAU2 = A(2)**2                                                    
!     TAU3 = A(3)**2                                                    
!  610 FORMAT(1H ,10F10.5 )                                              
!                                                                       
para(4) = tau1
para(5) = tau2
para(6) = tau3
if( m2 .eq. 0 )  go to 40
do 30  i=1,m2
!xx   30 atmp(I) = 0.90D0*DSIN( A(L+I) )                                   
atmp(i) = 0.90d0*dsin( a(l+i) )
30 continue
!c      CALL  ARCOEF( atmp,M2,A2 )
call  arcoefd( atmp,m2,ar )
do 35 i=1,m2
!c 35      para(i+6) = A2(i)
!xx 35      para(i+6) = AR(i)
para(i+6) = ar(i)
35 continue
40 continue
do 41 i=1,7
!xx 41      para(i+6+m2) = TDF(i)
para(i+6+m2) = tdf(i)
41 continue
!                                                                       
return
end
