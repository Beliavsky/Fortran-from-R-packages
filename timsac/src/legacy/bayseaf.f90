! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

! Reviced	M.S C85-02-19-16:46:06 BAYSEAA PAIR                         
!C      PROGRAM BAYSEA
subroutine bayseaf(y,ndata,focast,cdata,dmoi,trend,season,tdcmp,&
&irreg,adjust,est,psds,psdt,avabic,ipara,para,&
&arft,arfs,arfn,iart,iars,iarn)
  use timsac_kinds, only: dp
  implicit none
!
!     ---      --       -                                               
!     BAYESIAN SEASONAL ADJUSTMENT PROCEDURE                            
!                                                                       
!     THIS IS VERSION(3/1/85) OF BAYSEA WHICH WAS ORIGINALLY            
!     PUBLISHED IN                                                      
!                                                                       
!       AKAIKE,H. AND ISHIGURO,M. (1980)                                
!         BAYSEA, A BAYESIAN SEASONAL ADJUSTMENT PROGRAM.               
!         COMPUTER SCIENCE MONOGRAPHS, NO.13,                           
!         THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO.              
!                                                                       
!     THIS VERSION WAS DESIGNED AND PROGRAMMED BY                       
!     HIROTUGU AKAIKE AND MAKIO ISHIGURO, THE INSTITUTE OF STATISTICAL  
!     MATHEMATICS, 4-6-7 MINAMI-AZABU, MINATO-KU, TOKYO, 106, JAPAN.    
!     THE SUBROUTINES FOR OUTLIER CORRECTION WERE PREPARED BY           
!     GENSHIRO KITAGAWA.                                                
!                                                                       
!     THIS PROGRAM REALIZES A DECOMPOSITION OF TIME SERIES Y            
!     INTO THE FORM                                                     
!     Y(I) = T(I) +S(I)+I(I)+TDC(I)+OCF(I)                              
!     WHERE  T(I)=TREND  S(I)=SEASONAL  I(I)=IRREGULAR                  
!            TDC(I)=TRADING DAY COMPONENT     AND                       
!            OCF(I)=OUTLIER CORRECTION FACTOR                           
!                                                                       
!     THE PROCEDURE IS BASED ON A BAYESIAN MODEL AND ITS                
!     PERFORMANCE IS CONTROLLED BY THE SELECTION OF THE PARAMETERS OF   
!     THE PRIOR DISTRIBUTION.  THE CONSTRUCTION OF THE BASIC MODEL IS   
!     DISCUSSED IN THE FOLLOWING PAPERS:                                
!                                                                       
!       AKAIKE, H. (1980) LIKELIHOOD AND THE BAYES PROCEDURE.           
!         BAYESIAN STATISTICS, J.M.BERNARDO, M.H.DE GROOT, D.V.LINDLEY  
!         AND A.F.M.SMITH, EDS., UNIVERSITY PRESS, VALENCIA,            
!         SPAIN, 143-166.                                               
!       AKAIKE, H. (1980) SEASONAL ADJUSTMENT BY A BAYESIAN MODELING.   
!         JOURNAL OF TIME SERIES ANALYSIS, 1, 1-13.                     
!       AKAIKE, H. AND ISHIGURO, M. (1980) TREND ESTIMATION WITH        
!         MISSING OBSERVATIONS. ANNALS OF THE INSTITUTE OF STATISTICAL  
!         MATHEMATICS, 32,B, 481-488.                                   
!       ISHIGURO, M. AND AKAIKE, H. (1980) A BAYESIAN APPROACH TO       
!         THE TRADING-DAY ADJUSTMENT OF MONTHLY DATA.                   
!         TIME SERIES ANALYSIS, O.D.ANDERSON AND M.R.PERRYMAN, EDS.,    
!         NORTH-HOLLAND, AMSTERDAM, 213-226.                            
!       ISHIGURO,M. (1984) COMPUTATIONALLY EFFICIENT IMPLEMENTATION     
!         OF A BYAESIAN SEASONAL ADJUSTMENT PROCEDURE.                  
!         JOURNAL OF TIME SERIES ANALYSIS (TO APPEAR).                  
!       KITAGAWA,G. AND AKAIKE,H. (1982) A QUASI BAYESIAN APPROACH      
!         TO OUTLIER DETECTION. ANNALS OF THE INSTITUTE OF              
!         STATISTICAL MATHEMATICS, 34,B, 389-398.                       
!                                                                       
!     THE PRIOR DISTRIBUTION CONTROLS THE SMOOTHNESS OF THE TREND AND   
!     SEASONAL COMPONENTS BY ASSUMING LOW ORDER GAUSSIAN AR-MODELS      
!     FOR SOME DIFFERENCES OF THESE COMPONENTS. THE CHOICE OF           
!     THE VARIANCE OF THE GAUSSIAN DISTRIBUTION IS                      
!     REALIZED BY MAXIMIZING THE LOG LIKELIHOOD OF THE BAYESIAN MODEL.  
!                                                                       
!     FOR THE PURPOSE OF COMPARISON OF MODELS WITH DIFFERENT STRUCTURES 
!     THE CRITERION ABIC IS DEFINED BY                                  
!                                                                       
!          ABIC = (-2)LOG MAXIMUM LIKELIHOOD OF THE MODEL.              
!                                                                       
!     SMALLER VALUE OF ABIC REPRESENTS BETTER FIT.                      
!     FOR THE COMPARISON OF OVERALL PERFORMANCES OF VARIOUS             
!     MODELS THE AVERAGED ABIC (AVABIC) IS USEFUL.                      
!                                                                       
!     THIS PROGRAM REQUIRES THE FOLLOWING PARAMETERS.  WITHOUT FURTHER  
!     SPECIFICATION BY THE PROGRAM USER, THEY ARE RESPECTIVELY SET EQUAL
!     TO THE VALUES GIVEN IN THE PARENTHESES AT THE ENDS OF THEIR       
!     DESCRIPTIONS.                                                     
!                                                                       
!     PARAMETERS:                                                       
!       LOGT:   LOG-ADDITIVE-MODEL OPTION (0)                           
!        IF( LOGT.EQ.0) ADDITIVE MODEL                                  
!        IF( LOGT.EQ.1) LOG-ADDITIVE-MODEL                              
!       MT:     INPUT DEVICE SPECIFICATION (5)                          
!       RLIM:   OUTLIER LIMIT (0.0)                                     
!        IF( RLIM.GT.0.0) ANY DATA WHOSE VALUE IS GREATER               
!                         THAN RLIM IS TREATED AS A MISSING OBSERVATION 
!        IF( RLIM.LE.0.0) NO MISSING OBSERVATIONS                       
!       PERIOD: NUMBER OF SEASONALS WITHIN A PERIOD (12)                
!       SPAN:   NUMBER OF PERIODS TO BE PROCESSED AT ONE TIME (4)       
!       SHIFT:  NUMBER OF PERIODS TO BE SHIFTED                         
!               TO DEFINE THE NEW SPAN OF DATA (1)                      
!       FOCAST: LENGTH OF FORECAST AT THE END OF DATA (0)               
!       ORDER:  ORDER OF DIFFERENCING OF TREND (2)                      
!       ARFT(I): I-TH PARCOR OF DIFFERENCED TREND (I LESS THAN 4)       
!       (0.D0)                                                          
!       SORDER: ORDER OF DIFFERENCING OF SEASONAL (1)                   
!       ARFS(I): I*PERIOD-TH PARCOR OF DIFFERENCD SEASONAL              
!       (I LESS THAN 4) (0.D0)                                          
!       ARFN(I): I-TH PARCOR OF DIFFERENCD SEASONAL                     
!       (I LESS THAN 4) (0.DO)                                          
!       (FOR THE INITIAL CHOICES OF ARFT,ARFS,AND ARFN, USE             
!        PARCOR (PARTIAL AUTOCORRELATION) OUTPUTS FOR THE               
!        DEFAULT OPTION)                                                
!               (SORDER SHOULD NOT BE GREATER THAN SPAN.)               
!       RIGID:  CONTROLS THE RIGIDITY OF THE SEASONAL COMPONENT, MORE   
!               RIGID SEASONAL WITH LARGER RIGID (1.0)                  
!               (TO BE ADJUSTED ONLY AFTER THE SELECTION OF ORDER AND   
!                SORDER)                                                
!       YEAR:   TRADING-DAY ADJUSTMENT OPTION  (0)                      
!        IF(YEAR .EQ. 0) WITHOUT TRADING-DAY ADJUSTMENT                 
!        IF(YEAR .NE. 0) WITH TRADING-DAY ADJUSTMENT                    
!         NOTE: THE SERIES IS SUPPOSED TO START AT THIS 'YEAR'          
!       MONTH:  NUMBER OF THE MONTH IN WHICH THE SERIES STARTS (1)      
!        IF(YEAR .EQ. 0) THIS PARAMETER IS IGNORED                      
!       WTRD:   CONTROLS THE ADAPTIVITY OF THE TRADING-DAY-COMPONENT    
!               MORE DATA-ADAPTIVE WITH SMALLER WTRD (1.0)              
!       IOUTD: OUTLIER CORRECTION OPTION (0)                            
!        IF(IOUTD .EQ. 0) WITHOUT OUTLIER DETECTION                     
!        IF(IOUTD .EQ. 1) WITH OUTLIER DETECTION BY MARGINAL PROBABILITY
!        IF(IOUTD .EQ. 2) WITH OUTLIER DETECTION BY MODEL SELECTION     
!         (NOTE: OUTLIER DETECTION IS                                   
!                EXPENSIVE. DECLARING ABNORMAL OBSERVATIONS AS MISSING  
!                BY RLIM IS MUCH CHEAPER)                               
!       SPEC:  SPECTRUM ESTIMATION OPTION (1)                           
!        IF(SPEC .EQ. 0) NO SPECTRUM                                    
!        IF(SPEC .EQ. 1) SPECTRA OF IRREGULAR AND                       
!                         DIFFERENCED ADJUSTED                          
!         (IF THERE ARE MISSING OBSERVATIONS, SPECTRA OF THE            
!          INTERPOLATED DATA ARE GIVEN)                                 
!       PUNCH:  CARD OUTPUT CONTROL (0)                                 
!        IF(PUNCH .EQ. 0) OUTPUT SUPRESSED                              
!        IF(PUNCH .EQ. 1) CARD OUTPUT AVAILABLE                         
!                                                                       
!      --- PARAMETERS BELOW THIS LINE ARE SELDOM TO BE MODIFIED ---     
!                                                                       
!       ZERSUM: CONTROLS THE SUM OF THE SEASONALS WITHIN A PERIOD,      
!               CLOSER TO ZERO WITH LARGER ZERSUM (1.0)                 
!       DELTA:  CONTROLS THE LEAP YEAR EFFECT (7.0)                     
!       ALPHA:  CONTROLS PRIOR VARIANCE OF INITIAL TREND (0.01D0)       
!       BETA:   CONTROLS PRIOR VARIANCE OF INITIAL SEASONAL             
!               (0.01D0)                                                
!       GAMMA:  CONTROLS PRIOR VARIANCE OF INITIAL SUM OF SEASONAL      
!               (0.01D0)                                                
!                                                                       
!     THE FOLLOWING PROCEDURE OF PARAMETER MODIFICATION IS RECOMMENDED  
!     FOR USUAL APPLICATIONS:                                           
!     FIRST TRY THE COMBINATION(ORDER=2,SORDER=1). IF NECESSARY, TRY    
!     (ORDER=2,SORDER=2). THEN TRY A REDUCED VALUE OF 'RIGID' AND CHECK 
!     ABIC. WHEN 'RIGID' IS SUITABLY CHOSEN, USUALLY SHALLOW TROUGHS    
!     APPEAR AT FRIQUENCIES MARKED BY '---X' ON THE CHARTS OF SPECTRA.  
!                                                                       
!     TO MODIFY THE PARAMETERS FOLLOW THE FOLLOWING EXAMPLE:            
!       IF THE USER WANTS TO SET SORDER=2,RIGID=0.5 AND SO ON, AND IF   
!       MT IS CARD READER, PLACE THE CARDS WHICH CONTAIN THE FOLLOWING  
!       THREE TYPES OF STATEMENTS ON TOP OF THE INPUT DATA( PUNCH ONE   
!       SPACE AT THE FIRST COLUMN OF EACH CARD) :                       
!  &PARAM  (ON THE FIRST CARD)                                          
!  SORDER=2, RIGID=0.5, AND SO ON, (ON THE SECOND AND LATER CARDS)      
!  &END  (ON THE LAST CARD)                                             
!     NOTE THE COMMA AT THE END OF EACH PARAMETER SPECIFICATION.        
!                                                                       
!     INPUT DATA:                                                       
!                                                                       
!     THE FOLLOWING DATA SET SHOULD BE FED THROUGH THE INPUT DEVICE     
!     SPECIFIED BY MT IN THE FORMATS SHOWN IN THE PARENTHSES.           
!                                                                       
!     TITLE (20A4):   TITLE OF THE DATA                                 
!     NDATA (I5):     DATA LENGTH                                       
!     FORMAT SPECIFICATION OF DATA (20A4): FOR EXAMPLE, (4D20.10)       
!     DATA:       Y(I);I=1,NDATA                                        
!                                                                       
!     NUMERICAL OUTPUTS:                                                
!       OCF:    OUTLIER CORRECTION FACTOR                               
!       TREND:  ESTIMATED TREND                                         
!       SEASONAL: ESTIMATED SEASONAL                                    
!       TDCMP: ESTIMATED TRADING-DAY COMPONENT                          
!       IRREGULAR = 0.0        (IF OBSERVATION IS MISSING)              
!                 = ORIGINAL DATA - TREND - SEASONAL - TDCMP - OCF      
!                              (OTHERWISE)                              
!       ADJUSTED = TREND + IRREGULAR                                    
!       SMOOTHED = TREND + SEASONAL + TDCMP                             
!                                                                       
!     GRAPHICAL OUTPUTS:                                                
!       ORIGINAL DATA                                                   
!       OCF                                                             
!       TREND                                                           
!       SEASONAL                                                        
!       IRREGULAR                                                       
!       ADJUSTED                                                        
!       SMOOTHED                                                        
!       TDCMP                                                           
!       MISSING OBSERVATION INTERPOLATED DATA                           
!       SPECTRA OF IRREGULAR AND DIFFERENCED ADJUSTED                   
!                                                                       
!     WORKING AREA REQUIRED BY THIS PROGRAM:                            
!     VECTORS:                                                          
!       Y(IY)                                                           
!       TREND(IRSLT)                                                    
!       SEASON(IRSLT)                                                   
!       EST(IRSLT)                                                      
!       IRREG(IRSLT)                                                    
!       TDCMP(IRSLT)                                                    
!       ADJUST(IRSLT)                                                   
!       CDATA(IRSLT)                                                    
!       DMOI(IRSLT)                                                     
!       FTRN(IA)                                                        
!       FSEA(IA)                                                        
!       YS(IA)                                                          
!       YS1(IA)                                                         
!       YO(IA)                                                          
!       DC(MDC)                                                         
!       H2(8*IA)                                                        
!       WEEK(7*IA)                                                      
!   WHERE                                                               
!    IY = MAXIMUM DATA LENGTH                                           
!    IRSLT = MAXIMUM OUTPUT LENGTH                                      
!    IA = ( MAXIMUM DATA LENGTH WITHIN A SPAN)                          
!                                                                       
!   STRUCTURE OF THE PROGRAM                                            
!                                                                       
!   ****************00002150
!   *              *00002160
!   *    BAYSEA    *00002170
!   *              *00002180
!   ********I*******00002190
!           I                                                           
!           I-ARCOEF---PARTAR                                           
!           I                                                           
!           I----------------------------DSQRT                          
!           I                                                           
!           I-DATAR ------------DLOG                                    
!           I                                                           
!           I-------------------DLOG                                    
!           I                                                           
!           I-------------------COPY                                    
!           I                                                           
!           I-CALEND------------MOD                                     
!           I                                                           
!           I-SUBSEA-+-------------------DSQRT                          
!           I        I                                                  
!           I        I-SETDC -+-SETD  ---MIN0                           
!           I        I        I                                         
!           I        I        I-INIT                                    
!           I        I        I                                         
!           I        I        +-EXHSLD                                  
!           I        I                                                  
!           I        I-------------------DABS                           
!           I        I                                                  
!           I        I----------DLOG                                    
!           I        I                                                  
!           I        I-SETX  -+-HUSHLD-+-DABS                           
!           I        I        I        I                                
!           I        I        I        +-DSQRT                          
!           I        I        I                                         
!           I        I        I----------DABS                           
!           I        I        I                                         
!           I        I        I-DLOG                                    
!           I        I        I                                         
!           I        I        I-SETD  ---MIN0                           
!           I        I        I                                         
!           I        I        I-INIT                                    
!           I        I        I                                         
!           I        I        I-EXHSLD                                  
!           I        I        I                                         
!           I        I        +-MOD                                     
!           I        I                                                  
!           I        I-SOLVE                                            
!           I        I                                                  
!           I        I-DECODE-+-CLEAR                                   
!           I        I        I                                         
!           I        I        I----------DSQRT                          
!           I        I        I                                         
!           I        I        I-COPY                                    
!           I        I        I                                         
!           I        I        I-PRDCT                                   
!           I        I        I                                         
!           I        I        I-ADD                                     
!           I        I        I                                         
!           I        I        +-SBTRCT                                  
!           I        I                                                  
!           I        +-DEXP                                             
!           I                                                           
!           I-------------------SBTRCT                                  
!           I                                                           
!           I-OUTLIR-+-SRTMIN                                           
!           I        I                                                  
!           I        I----------DLOG                                    
!           I        I                                                  
!           I        I----------DFLOAT                                  
!           I        I                                                  
!           I        I----------BINARY                                  
!           I        I                                                  
!           I        I-LKOUT1-+-DFLOAT                                  
!           I        I        I                                         
!           I        I        I-DLOG                                    
!           I        I        I                                         
!           I        I        I-POOLAV                                  
!           I        I        I                                         
!           I        I        I----------DSQRT                          
!           I        I        I                                         
!           I        I        +-PERMUT---ISORT                          
!           I        I                                                  
!           I        I-DEXP                                             
!           I        I                                                  
!           I        I-PRPOST---BINARY                                  
!           I        I                                                  
!           I        +-MODIFY---BINARY                                  
!           I                                                           
!           I-------------------ADD                                     
!           I                                                           
!           I----------DEXP                                             
!           I                                                           
!           I-DFR1                                                      
!           I                                                           
!           I-GRAPH -+-LOG10                                            
!           I        I                                                  
!           I        I-SQRT                                             
!           I        I                                                  
!           I        I-------------------ABS                            
!           I        I                                                  
!           I        I----------MOD                                     
!           I        I                                                  
!           I        +-DEXP                                             
!           I                                                           
!           +-SPGRH -+-------------------MIN0                           
!                    I                                                  
!                    I-SAUTCO-+-SMEADL---SUMF                           
!                    I        I                                         
!                    I        I-CROSCO---DBLE                           
!                    I        I                                         
!                    I        +-CORNOM---DSQRT                          
!                    I                                                  
!                    I-------------------DSQRT                          
!                    I                                                  
!                    I-SICP2 -+-DLOG                                    
!                    I        I                                         
!                    I        +----------DSQRT                          
!                    I                                                  
!                    +-SNRASP-+-FOUGER-+-DCOS                           
!                             I        I                                
!                             I        +-DSIN                           
!                             I                                         
!                             I-SUBVCP                                  
!                             I                                         
!                             I-DLOG10                                  
!                             I                                         
!                             I-DSP3  -+-AMAX                           
!                             I        I                                
!                             I        I-ABS                            
!                             I        I                                
!                             I        +-AMIN                           
!                             I                                         
!                             +-MOD                                     
!                                                                       
!c      !DEC$ ATTRIBUTES DLLEXPORT::BAYSEAF
!
!xx       IMPLICIT  REAL*8( A-H,O-Z )                                      
!xx      INTEGER*4  ORDER, SORDER, PERIOD, SPAN, OVLAP, FOCAST, HEAD, SHIFT
!c     * ,TAIL,PUNCH,YEAR,SPEC                                            
!xx     * ,TAIL,YEAR,SPEC                                            
!xx      REAL*8 IRREG, IRREG0                                              
!c      DIMENSION  FTRN(200),FSEA(200),PSDS(500),PSDT(500),PSDS0(500),    
!c     *           PSDT0(500)                                             
!c      DIMENSION  SEASON(500), TREND(500), EST(500), ADJUST(500),        
!c     * IRREG(500),WEEK(2000),TDCMP(500),TDCMP0(500),DMOI(500)           
!c     *       ,SEAS0(500), TREND0(500), EST0(500), ADJ0(500), IRREG0(500)
!c      DIMENSION  DC(40000),ARFT(3),ARFS(3),ARFN(3), H2(4000)            
!c      DIMENSION  Y(500),YS(500),CDATA(500),YS1(500),YO(500)             
!c      DIMENSION   DADJ(500),H(20000),F(1000)
!c      DIMENSION  DC(40000),H2(4000),H(20000),F(1000),WEEK(2000)
!c      DIMENSION  YS(500),YS1(500),YO(500)
!xx      DIMENSION  Y(NDATA),CDATA(NDATA),DMOI(NDATA),
!xx     *            TREND(NDATA+FOCAST),SEASON(NDATA+FOCAST),
!xx     *            TDCMP(NDATA+FOCAST),IRREG(NDATA),
!xx     *            ADJUST(NDATA),EST(NDATA+FOCAST),
!xx     *            PSDS(NDATA+FOCAST),PSDT(NDATA+FOCAST),
!xx     *            IPARA(12),PARA(8),ARFT(3),ARFS(3),ARFN(3)
!xx      DIMENSION  TREND0(NDATA+FOCAST),SEAS0(NDATA+FOCAST),
!xx     *            TDCMP0(NDATA+FOCAST),IRREG0(NDATA+FOCAST),
!xx     *            ADJ0(NDATA+FOCAST),EST0(NDATA+FOCAST),
!xx     *            PSDS0(NDATA+FOCAST),PSDT0(NDATA+FOCAST)
!xx      DIMENSION  FTRN(IPARA(4)+3),FSEA((IPARA(5)+3)*IPARA(1)+3)
!xx      DIMENSION  F(NDATA+FOCAST+1),WEEK(7,NDATA+FOCAST)
!xx      DIMENSION  YS(NDATA),YS1(NDATA),YO(NDATA)
!
integer ndata, focast, ipara(12), iart, iars, iarn
real(dp) y(ndata), cdata(ndata), dmoi(ndata),&
&trend(ndata+focast), season(ndata+focast),&
&tdcmp(ndata+focast), irreg(ndata), adjust(ndata),&
&est(ndata+focast), psds(ndata+focast),&
&psdt(ndata+focast), avabic, para(8), arft(3),&
&arfs(3), arfn(3)
! local
integer i, is, iq, iend, icnt, icnt1, iout, ioutd, idc, istem,&
&itrn, i1, i2, n, n2, npf, nh, nday, next, ntem,&
&l, logt, lf, lftrn, lfsea, limit, linkt, links, month, nf,&
&order, sorder, period, span, ovlap, head,&
&shift, tail, year, spec
real(dp) trend0(ndata+focast), seas0(ndata+focast),&
&tdcmp0(ndata+focast), irreg0(ndata+focast),&
&adj0(ndata+focast), est0(ndata+focast),&
&psds0(ndata+focast), psdt0(ndata+focast),&
&ftrn(ipara(4)+3), fsea((ipara(5)+3)*ipara(1)+3),&
&f(ndata+focast+1), week(7,ndata+focast),&
&ys(ndata), ys1(ndata), yo(ndata), rlim, rigid,&
&wtrd, dd, zersum, delta, alpha, beta, gamma,&
&an, ap, zer, smth, smth2, rout, sy, ytem, count,&
&abic
!c      CHARACTER*80   TITLE
!c      COMMON /ILOGT/ LOGT,ISHRNK,PUNCH,IOUTD,ROUT                       
!c      COMMON /IDATA/ PERIOD,ORDER,SORDER,YEAR,NDAY,IFIX                 
!c      COMMON /RDATA/ ALPHA,BETA,GAMMA,ZER,SMTH,SMTH2,DD,WTRD,DELTA      
!c      NAMELIST /PARAM/  MT, PERIOD, RLIM, SPAN, SHIFT, FOCAST,          
!c     *                  RIGID, ORDER, SORDER, ZERSUM,  LOGT, PUNCH      
!c     * ,YEAR,MONTH,WTRD,DELTA,SPEC,IOUTD,IFIX,ARFT,ARFS,ARFN,IART       
!c     * ,IARS,IARN
!
!      FILE NAME INPUT 
!
!c 	CHARACTER(100)  IFLNAM,OFLNAM,MFLNAM
!c         CALL FLNAM3(IFLNAM,OFLNAM,MFLNAM,NFL)
!c            NFL = 4
!c 	IF (NFL.EQ.0) GO TO 999
!c 	IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) THEN
!c 	   OPEN (6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR)
!c 	ELSE
!c 	   CALL SETWND
!c 	END IF
!                                                                       
npf = ndata+focast
!                                                                       
!     PRESET CONSTANTS                                                  
!                                                                       
!C 1496 REWIND 5                                                          
!c      IY = 500                                                          
!c      IA = 500                                                       
!c      IRSLT = 500                                                       
!c      MDC=40000                                                         
!c      ALPHA = 0.01D0                                                    
!c      BETA = 0.01D0                                                     
!c      GAMMA = 0.1D0                                                     
!                                                                       
!     PRESET PARAMETERS                                                 
!                                                                       
!c      RIGID = 1.D0                                                      
!c      ZERSUM = 1.D0                                                     
!c      ORDER = 2                                                         
!c      SORDER = 1                                                        
!c      PERIOD = 12                                                       
!c      SPAN = 4                                                          
!c      SHIFT = 1                                                         
!c      FOCAST = 0                                                        
rlim = 0.0
!c      MT = 5                                                            
!c      DD = 1.D0                                                         
!c      LOGT=0                                                            
!c      PUNCH = 0                                                         
!c      WTRD = 1.D0                                                       
!c      YEAR = 0                                                          
!c      MONTH = 1                                                         
!c      IOUTD=0                                                           
!c      SPEC = 1                                                          
!c      NDAY = 1                                                          
!c      DELTA=7.D0                                                        
period=ipara(1)
span=ipara(2)
shift=ipara(3)
order=ipara(4)
sorder=ipara(5)
logt=ipara(6)
year=ipara(7)
month=ipara(8)
nday=ipara(9)
spec=ipara(10)
ioutd=ipara(11)
idc=ipara(12)
rigid=para(1)
wtrd=para(2)
dd= para(3)
zersum=para(4)
delta=para(5)
alpha=para(6)
beta=para(7)
gamma=para(8)
!c      ARFT(1) = 0.0D0                                                   
!c      ARFT(2) = 0.D0                                                    
!c      ARFT(3) = 0.D0                                                    
!c      ARFS(1) = 0.0D0                                                   
!c      ARFS(2) = 0.0D0                                                   
!c      ARFS(3) = 0.D0                                                    
!c      ARFN(1) = 0.0D0                                                   
!c      ARFN(2) = 0.0D0                                                   
!c      ARFN(3) = 0.D0                                                    
!c      WRITE(6,1)                                                        
!c    1 FORMAT(1H ,'BAYSEA',/                                             
!c     * 1H ,'BAYESIAN TREND AND SEASONAL ESTIMATION OF TIME SERIES'      
!c     * ,'  VERSION 3/1/85',/                                            
!c     * ,'   WITH TRADING-DAY AND LEAP-YEAR ADJUSTMENT AND',/            
!c     * ,'        OUTLIER DETECTION OPTION')                             
!                                                                       
!     PARAMETER MODIFICATION AND                                        
!     DATA INPUT                                                        
!                                                                       
!c      READ(5,PARAM)
if(   sorder .gt. span ) sorder = span
!c      NDAY=1                                                            
!c      WRITE(6,610) LOGT,MT,RLIM,PERIOD,SPAN,SHIFT,FOCAST,ORDER,         
!c     *             SORDER,RIGID,YEAR,MONTH,WTRD,IOUTD,SPEC,PUNCH,       
!c     *             ZERSUM,DELTA
!c      CALL ARCOEF('ARFT',ARFT,IART)                                     
!c      CALL ARCOEF('ARFS',ARFS,IARS)                                     
!c      CALL ARCOEF('ARFN',ARFN,IARN)                       
is = period*sorder
ap=period
!c      IPRD=2                                                            
!c      IF(PERIOD.EQ.1)IPRD=1                                             
lftrn = order + iart
lfsea = (sorder + iars)*period + iarn
!c      IDC=LFTRN*IPRD+1                                                  
!c      IDCX = LFSEA * 2 + 1                                              
!c      IF(PERIOD .GT. 1 .AND. IDC .LT. IDCX)  IDC=IDCX                   
!c      IF(PERIOD .GT. 1 .AND. IDC .LT. PERIOD*2-1) IDC=PERIOD*2-1
nh= lfsea + 1
n2=1
if(year .ne. 0) n2=8
!  ************                                                         
zer=zersum/dsqrt(ap)*rigid
smth = 1.d0/rigid
smth2=1.d0
!  ************                                                         
!C      IF( MT .GT. 7 ) REWIND MT
!c      OPEN (MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD')
!c      IF (PUNCH.EQ.1) THEN
!c         IF ((NFL.EQ.3) .OR. (NFL.EQ.4)) THEN
!c            OPEN (7,FILE=MFLNAM,ERR=920,IOSTAT=IVAR)
!c         ELSE
!c            OPEN (7,FILE='baysea.out',ERR=930,IOSTAT=IVAR)
!c         END IF
!c      END IF
!
!c      CALL DATAR( TITLE,Y,IY,NDATA,MT,LOGT,RLIM )
!c      CLOSE(MT)
if(ioutd .eq. 0) go to 1212
rout = 1.d60
rlim = 1.d50
if(logt .eq. 0) go to 1212
rout = dlog(rout)
rlim = dlog(rlim)
!-----
do 20 i=1,ndata
!xx   20 Y(I) = DLOG(Y(I))
y(i) = dlog(y(i))
20 continue
!-----
1212 continue
!     -------------------                                               
!                                                                       
!                                                                       
!     WORKING AREA  CHECK                                               
!                                                                       
!c      NPAR=(SPAN*2-1)*PERIOD                                            
!c      IF(NPAR.GT.NDATA+FOCAST)NPAR=NDATA+FOCAST                         
!c      IF( NDATA+FOCAST .LE. IRSLT )   GO TO 100                         
!c      WRITE( 6,603 )                                                    
!c  603 FORMAT( 1H , 'IRSLT IS TOO SMALL OR NDATA+FOCAST IS TOO LARGE' )  
!c  100 IF( NPAR.LE.IA )   GO TO 200                                      
!c      WRITE( 6,605 )                                                    
!c  604 FORMAT( 1H ,'MDC IS TOO SMALL')                                   
!c  200 CONTINUE                                                          
!c      NPAR=NPAR*IPRD                                                    
!c      IF( NPAR .LE. MDC )   GO TO 300                                   
!c      WRITE( 6,604 )                                                    
!c  605 FORMAT( 1H ,'IA IS TOO SMALL' )                                   
!c  300 CONTINUE                                                          
!                                                                       
!                                                                       
!     INITIALIZATION                                                    
nf=order
if(period.gt.1.and.nf.lt.is)nf=is
!c      WRITE(6,606)                                                      
!c  606 FORMAT(1H ,'INITIALIZATION'    )                                  
!     N=LENGTH OF A SPAN                                                
ovlap = span - 1
limit = ndata - ovlap*period
sy = 0.d0
!-----
iq = 0
!-----
do 2468 i=1,ndata
ytem = y(i)
if(rlim .le. 0.d0) go to 4681
if(ytem .ge. rlim) go to 2468
4681 iq = iq + 1
sy = sy + ytem
if(iq .ge. period) go to 4680
2468 continue
4680 continue
ytem = sy / ap
do 8  i=1,lftrn
!xx    8 FTRN(I) = YTEM                                                    
ftrn(i) = ytem
8 continue
if(lfsea .eq. 0) go to 998
!xx      DO 9  I=1,LFSEA                                                   
!xx    9 FSEA(I) = 0.D0
fsea(1:lfsea) = 0.d0
998 continue
n = (span*2-1)*period
avabic = 0.d0
count = 0.d0
iend = 0
!                                                                       
!     ***************                                                   
!     **           **                                                   
!     ** MAIN LOOP **                                                   
!     **           **                                                   
!     ***************                                                   
!     ICNT0: ITERATION CONTROL FOR THE 0-TH SPAN                        
!                                                                       
do 1000  icnt1=1,1000
!                                                                       
!     MANIPULATION OF THE ICNT-TH SPAN                                  
!                                                                       
icnt = icnt1-1
!                                                                       
!     DATA END DETECTION                                                
!                                                                       
!     --------------------                                              
head = 1 + (icnt1-2)*shift*period+span*period
if(icnt1 .eq. 1) head=1
if( icnt .le. 0 )   go to 2345
!     --------------------                                              
!     HEAD: INITIAL POINT OF THE NEW SPAN                               
if(head.le.limit)   go to 2345
go to 1234
!     NEW SPAN                                                          
2345 continue
!     TAIL: END POINT OF THE NEW SPAN                                   
tail = head+n-1
if( tail .gt. ndata )   n = ndata-head+1
tail = head + n-1
!     --------------------                                              
if( tail .eq. ndata )   iend = 1
!     IEND=1: LAST SPAN IN NORMAL ITERATION                             
!      --------------------                                             
!                                                                       
!                                                                       
!c      CALL COPY(YO     ,N,1,N,1,1,Y     ,N,1,NDATA,HEAD,1)           
!xx         CALL BCOPY(YO     ,N,1,N,1,1,Y     ,N,1,NDATA,HEAD,1)
call bcopy(yo     ,n,1,1,1,y     ,n,1,head,1)
!                                                                       
!     INITIALIZATION FOR THE INNER LOOP SUBSEA                          
!                                                                       
!xx 6789 CONTINUE                                                          
if(year .ne. 0) call calend(week,year,month+head-1,n+focast)
!c      WRITE(6,2) ICNT1,HEAD,TAIL                                        
!xx    2 FORMAT(1H ,'(',I5,' )TH SPAN  HEAD =',I6,'   TAIL =',I6 )         
!c      WRITE( 6,620 )   (FTRN(I),I=1,LFTRN)                              
!xx  620 FORMAT(1H ,'* INITIAL TREND  *'/,(1X,12D11.3))                    
!c      IF(LFSEA .NE. 0) WRITE( 6,621 )   (FSEA(I),I=1,LFSEA)             
!xx  621 FORMAT(1H ,'* INITIAL SEASONAL  *'/,(1X,12D11.3))                 
next = head + n
if( iend .eq. 0 ) next = next - ovlap*period
linkt = next - lftrn
links = next - lfsea
!     --------------------                                              
!     DO SEARCH CONTROL                                                 
itrn = icnt1
!     --------------------                                              
iout=0
!c      CALL COPY(YS,N,1,N,1,1,YO, N,1,N,1,1)                             
!xx      CALL BCOPY(YS,N,1,N,1,1,YO, N,1,N,1,1)
call bcopy(ys,n,1,1,1,yo, n,1,1,1)
9700 continue
!     SEASONAL DECOMPOSITION OF ICN-TH LOCAL SPAN                       
call subsea(abic,seas0,trend0,est0,adj0,irreg0,tdcmp0,&
!c     *  FSEA,FTRN,YS,N,FOCAST,RLIM,WEEK,DC,IDC,H,NH,F,H2,N2,ITRN,       
!c     *  IARS,ARFS,IART,ARFT,IARN,ARFN,PSDT0,PSDS0)                      
&fsea,lfsea,ftrn,ys,n,focast,rlim,week,idc,nh,f,n2,itrn,&
&iars,arfs,iart,arft,iarn,arfn,psdt0,psds0,npf,&
&period,order,sorder,year,nday,logt,&
&alpha,beta,gamma,zer,smth,smth2,dd,wtrd,delta)
if(ioutd .eq. 0) go to 9600
if(iout .ge. 2) go to 9600
call sbtrct(irreg0,n,yo,n,est0,n)
!c      CALL OUTLIR(IRREG0,N,10,2,1,YS1,RLIM) 
call outlir(irreg0,n,10,2,1,ys1,rlim,ioutd,rout)
iout=iout+1
call add(ys,n,est0,n,ys1,n)
!c      IF(IOUT .EQ. 2) CALL COPY(YS1,N,1,N,1,1, YS,N,1,N,1,1)     
!xx      IF(IOUT .EQ. 2) CALL BCOPY(YS1,N,1,N,1,1, YS,N,1,N,1,1)
if(iout .eq. 2) call bcopy(ys1,n,1,1,1, ys,n,1,1,1)
if(iout .ne. 1) go to 9700
ntem = n+1
do 9702 i=1,order
ntem=ntem-1
!xx 9702 YS(NTEM)=ROUT                                                     
ys(ntem)=rout
9702 continue
go to 9700
9600 continue
!                                                                       
!                                                                       
!     RECORDING THE BEST RESULT                                         
!                                                                       
l=next - head
lf=l+focast
!c      CALL COPY(PSDT,LF,1,LF,HEAD,1,PSDT0,LF,1,LF,1,1)                  
!c      CALL COPY(PSDS,LF,1,LF,HEAD,1,PSDS0,LF,1,LF,1,1)                  
!c      CALL COPY(SEASON,LF,1,LF,HEAD,1,SEAS0,LF,1,LF,1,1)                
!c      CALL COPY(TREND,LF,1,LF,HEAD,1,TREND0,LF,1,LF,1,1)                
!c      CALL COPY(EST,LF,1,LF,HEAD,1,EST0,LF,1,LF,1,1)                    
!c      CALL COPY(ADJUST,L,1,L,HEAD,1,ADJ0,L,1,L,1,1)                     
!c      CALL COPY(IRREG,L,1,L,HEAD,1,IRREG0,L,1,L,1,1)                    
!c      CALL COPY(TDCMP,LF,1,LF,HEAD,1,TDCMP0,LF,1,LF,1,1)                
!c      CALL COPY(CDATA,L,1,L,HEAD,1,YS1,L,1,L,1,1)         
!xx      CALL BCOPY(PSDT,LF,1,LF,HEAD,1,PSDT0,LF,1,LF,1,1)
!xx      CALL BCOPY(PSDS,LF,1,LF,HEAD,1,PSDS0,LF,1,LF,1,1)
!xx      CALL BCOPY(SEASON,LF,1,LF,HEAD,1,SEAS0,LF,1,LF,1,1)
!xx      CALL BCOPY(TREND,LF,1,LF,HEAD,1,TREND0,LF,1,LF,1,1)
!xx      CALL BCOPY(EST,LF,1,LF,HEAD,1,EST0,LF,1,LF,1,1)
!xx      CALL BCOPY(ADJUST,L,1,L,HEAD,1,ADJ0,L,1,L,1,1)
!xx      CALL BCOPY(IRREG,L,1,L,HEAD,1,IRREG0,L,1,L,1,1)
!xx      CALL BCOPY(TDCMP,LF,1,LF,HEAD,1,TDCMP0,LF,1,LF,1,1)
!xx      CALL BCOPY(CDATA,L,1,L,HEAD,1,YS1,L,1,L,1,1)
call bcopy(psdt,lf,1,head,1,psdt0,lf,1,1,1)
call bcopy(psds,lf,1,head,1,psds0,lf,1,1,1)
call bcopy(season,lf,1,head,1,seas0,lf,1,1,1)
call bcopy(trend,lf,1,head,1,trend0,lf,1,1,1)
call bcopy(est,lf,1,head,1,est0,lf,1,1,1)
call bcopy(adjust,l,1,head,1,adj0,l,1,1,1)
call bcopy(irreg,l,1,head,1,irreg0,l,1,1,1)
call bcopy(tdcmp,lf,1,head,1,tdcmp0,lf,1,1,1)
call bcopy(cdata,l,1,head,1,ys1,l,1,1,1)
an = n
avabic = avabic + abic
count = count + an
!                                                                       
!                                                                       
!     INITIAL VALUES FOR THE NEXT SPAN                                  
!                                                                       
if(iend .eq. 1) go to 1234
!c      CALL  COPY( FTRN,LFTRN,1,LFTRN,1,1,TREND,LFTRN,1,IOUT,            
!xx      CALL  BCOPY( FTRN,LFTRN,1,LFTRN,1,1,TREND,LFTRN,1,IOUT,
!xx     * LINKT,1)                                                         
call  bcopy( ftrn,lftrn,1,1,1,trend,lftrn,1,linkt,1)
istem = lfsea
if(links .ge. 1) go to 1111
links=1-links
istem=istem-links
i1=lfsea+1
do 2222 i=1,links
i1=i1-1
i2=i1-istem
!xx 2222 FSEA(I1)=FSEA(I2)                                                 
fsea(i1)=fsea(i2)
2222 continue
links=1
1111 continue
!c      CALL COPY(FSEA,ISTEM,1,ISTEM,1,1,SEASON,ISTEM,1,IOUT,LINKS,1)     
!xx      CALL BCOPY(FSEA,ISTEM,1,ISTEM,1,1,SEASON,ISTEM,1,IOUT,LINKS,1)
call bcopy(fsea,istem,1,1,1,season,istem,1,links,1)
if(icnt1 .gt. 1) go to 1000
alpha = 1.d0
beta = 1.d0
gamma = 1.d0
n = span*period
if(n .gt. ndata) n=ndata
ovlap = span-shift
limit = ndata-ovlap*period
1000 continue
!     ************************                                          
!     *                      *                                          
!     * END OF THE MAIN LOOP *                                          
!     *                      *                                          
!     ************************                                          
!                                                                       
1234 continue
!                                                                       
!     NUMERICAL OUTPUTS                                                 
!                                                                       
do 4444 i=1,ndata
cdata(i)=y(i) - cdata(i)
dmoi(i) = y(i)
if(rlim .le. 0.d0) go to 4444
if(y(i) .lt. rlim.and.ioutd.eq.0) go to 4444
if(ioutd .ne. 0 .and. y(i) .gt. rlim) go to 4442
if(ioutd.ne.0.and.-cdata(i).lt.rlim)go to 4443
4442 continue
adjust(i)=trend(i)
dmoi(i)=est(i)
irreg(i)=0.d0
4443 cdata(i) = y(i) - dmoi(i)
4444 continue
!c      NPF = NDATA+FOCAST                                                
if(logt .eq. 0) go to 1250
do 1240 i=1,npf
trend(i) = dexp(trend(i))
season(i) = dexp(season(i))
est(i) = dexp(est(i))
tdcmp(i) = dexp(tdcmp(i))
if(i .gt. ndata) go to 1240
irreg(i) = dexp(irreg(i))
y(i) = dexp(y(i))
adjust(i) = dexp(adjust(i))
cdata(i) = dexp(cdata(i))
dmoi(i) = dexp(dmoi(i))
1240 continue
if(rlim .gt. 0.d0) rlim = dexp(rlim)
1250 continue
!c      DO 1324 I=1,NDATA                                                 
!c 1324 DADJ(I)=ADJUST(I)                                                 
!c      CALL DFR1(1,1,NDATA,DADJ,NDATA1)                                  
!c      IF(PUNCH .EQ. 0) GO TO 5555                                       
!c      WRITE(7,714) (TITLE(I),I=1,20)                                    
!c      WRITE(7,720) NDATA                                                
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (Y(I),I=1,NDATA)                                     
!c      WRITE(7,716)                                                      
!c      WRITE(7,720) NPF                                                  
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (TREND(I),I=1,NPF)                                   
!c      WRITE(7,717)                                                      
!c      WRITE(7,720) NPF                                                  
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (SEASON(I),I=1,NPF)                                  
!c      WRITE(7,718)                                                      
!c      WRITE(7,720) NDATA                                                
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (IRREG(I),I=1,NDATA)                                 
!c      WRITE(7,719)                                                      
!c      WRITE(7,720) NDATA                                                
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (ADJUST(I),I=1,NDATA)                                
!c      WRITE(7,703)                                                      
!c      WRITE(7,720) NDATA1                                               
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (DADJ(I),I=1,NDATA1)                                 
!c      WRITE(7,715)                                                      
!c      WRITE(7,720) NPF                                                  
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (EST(I),I=1,NPF)                                     
!c      IF(IOUTD .EQ. 0) GO TO 8500                                       
!c      WRITE(7,721)                                                      
!c      WRITE(7,720)NDATA                                                 
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (CDATA(I),I=1,NDATA)                                 
!c 8500 CONTINUE                                                          
!c      IF(YEAR .EQ. 0) GO TO 5552                                        
!c      WRITE(7,730)                                                      
!c      WRITE(7,720) NPF                                                  
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (TDCMP(I),I=1,NPF)                                   
!c 5552 IF(RLIM .LE. 0.D0) GO TO 5555                                     
!c      WRITE(7,777)                                                      
!c      WRITE(7,720) NDATA                                                
!c      WRITE(7,700)                                                      
!c      WRITE(7,702) (DMOI(I),I=1,NDATA)                                  
!c 5555 CONTINUE                                                          
avabic = avabic/count
avabic = avabic*ndata
!c      WRITE(6,600) AVABIC                                               
!c      IF(IOUTD .EQ. 0) GO TO 5553                                       
!c      WRITE(6,641)                                                      
!c      WRITE(6,602) (CDATA(I),I=1,NDATA)                                 
!c 5553 CONTINUE                                                          
!c      WRITE( 6,616 )                                                    
!c      WRITE( 6,602 )   (TREND(I),I=1,NPF)                               
!c      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
!c      WRITE( 6,617 )                                                    
!c      WRITE( 6,602 )   (SEASON(I),I=1,NPF)                              
!c      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
!c      IF(YEAR .EQ. 0) GO TO 5432                                        
!c      WRITE(6,640)                                                      
!c      WRITE(6,602) (TDCMP(I),I=1,NPF)                                   
!c      IF(FOCAST .NE. 0) WRITE(6,630) FOCAST                             
!c 5432 CONTINUE                                                          
!c      WRITE( 6,618 )                                                    
!c      WRITE( 6,602 )   (IRREG(I),I=1,NDATA)                             
!c      WRITE( 6,619 )                                                    
!c      WRITE( 6,602 )   (ADJUST(I),I=1,NDATA)                            
!c      WRITE( 6,615 )                                                    
!c      WRITE( 6,602 )   (EST(I),I=1,NPF)                                 
!c      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
!                                                                       
!     GRAPHICAL OUTPUTS                                                 
!                                                                       
!c      ISD = 1                                                           
!c      IF(LOGT .NE. 0) ISD = 2                                           
!c      WRITE( 6,611 )                                                    
!c      WRITE( 6,614 )   (TITLE(I),I=1,20)                                
!c      CALL GRAPH(0,PSDT,Y,NDATA,RLIM,1)                                 
!c      IF(RLIM .LE. 0.D0) GO TO 8301                                     
!c      WRITE(6,611)                                                      
!c      WRITE(6,776)                                                      
!c      CALL GRAPH(0,PSDT,DMOI,NDATA,RLIM,0)                              
!c 8301 CONTINUE                                                          
!c      WRITE( 6,611 )                                                    
!c      WRITE( 6,616 )                                                    
!c      CALL GRAPH(ISD,PSDT,TREND,NDATA+FOCAST,RLIM,0)                    
!c      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
!c      WRITE( 6,611 )                                                    
!c      WRITE( 6,619 )                                                    
!c      CALL GRAPH(0,PSDT,ADJUST,NDATA,RLIM,0)                            
!c      WRITE( 6,611 )                                                    
!c      WRITE( 6,615 )                                                    
!c      CALL GRAPH(0,PSDT,EST,NDATA+FOCAST,RLIM,0)                        
!c      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
!c      ITEM = -1                                                         
!c      IF(LOGT .NE. 0) ITEM=1                                            
!c      WRITE( 6,611 )                                                    
!c      WRITE( 6,617 )                                                    
!c      CALL GRAPH(ISD,PSDS,SEASON,NDATA+FOCAST,RLIM,ITEM)                
!c      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
!c      WRITE( 6,611 )                                                    
!c      WRITE( 6,618 )                                                    
!c      CALL GRAPH(0,PSDT,IRREG,NDATA,RLIM,0)                             
!c      IF(IOUTD .EQ. 0) GO TO 8300                                       
!c      WRITE(6,611)                                                      
!c      WRITE(6,641)                                                      
!c      CALL GRAPH(0,PSDT,CDATA,NDATA,RLIM,0)                             
!c 8300 CONTINUE                                                          
!c      IF(LOGT .EQ. 1) GO TO 3332                                        
!c      WRITE( 6,611 )                                                    
!c      WRITE( 6,622 )                                                    
!c      CALL  GRAPH(0,PSDT, IRREG,NDATA,RLIM,2 )                          
!c 3332 IF(YEAR .EQ. 0) GO TO 3333                                        
!c      WRITE(6,611)                                                      
!c      WRITE(6,640)                                                      
!c      CALL GRAPH(0,PSDT,TDCMP,NDATA+FOCAST,RLIM,0)                      
!c      IF(FOCAST .NE. 0) WRITE(6,630) FOCAST                             
!c 3333 CONTINUE                                                          
!c      IF(SPEC .EQ. 0) GO TO 9999                                        
!c      IF(PUNCH .NE. 0) WRITE(7,8888)                                    
!c 8888 FORMAT('IRREGULAR(SPECTRUM)')                                     
!c      WRITE(6,800)                                                      
!c      CALL SPGRH(IRREG,NDATA,1)                                         
!c      WRITE(6,810)                                                      
!c      IF(PUNCH .NE. 0) WRITE(7,8887)                                    
!c 8887 FORMAT('DADJ(SPECTRUM)')                                          
!c      CALL SPGRH(DADJ,NDATA1,1)                                         
!c  800 FORMAT(1H ,'SPECTRUM OF IRREGULAR')                               
!c  810 FORMAT(1H ,'SPECTRUM OF DIFFERENCED ADJUSTED SERIES')             
!c 1810 FORMAT(1H ,'PARCOR OF',I2,' TIME(S) DIFFERENCED TREND SERIES')    
!c 1811 FORMAT(1H ,'PARCOR OF',I2,' TIME(S) DIFFERENCED SEASONAL SERIES') 
!c 9999 CONTINUE                                                          
!c      WRITE(6,1810) ORDER                                               
!c      CALL DFR1(ORDER,1,NDATA,TREND,NDATA1)                             
!c      CALL SPGRH(TREND,NDATA1,0)                                        
!c      WRITE(6,1811) SORDER                                              
!c      CALL DFR1(SORDER,PERIOD,NDATA,SEASON,NDATA1)                      
!c      CALL SPGRH(SEASON,NDATA1,0)
!c      GO TO 990
!
!c  900 CONTINUE
!c      WRITE(6,690) IVAR,OFLNAM
!c      GO TO 999
!c  910 CONTINUE
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
!c      WRITE(6,691) IVAR,IFLNAM
!c      GO TO 999
!c  920 CONTINUE
!c      CLOSE(5)
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
!c      WRITE(6,692) IVAR,MFLNAM
!c      GO TO 999
!c  930 CONTINUE
!c      CLOSE(5)
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
!c      WRITE(6,693) IVAR
!c      GO TO 999
!
!xx  690 FORMAT(1H ,' !!! Output_Data_File OPEN ERROR ',I8/1H ,100A)
!xx  691 FORMAT(1H ,' !!! Input_Data_File OPEN ERROR ',I8/1H ,100A)
!xx  692 FORMAT(1H ,' !!! Intermediate_Data_File OPEN ERROR ',I8/1H ,100A)
!xx  693 FORMAT(1H ,' !!! baysea.out  OPEN ERROR ',I8)
!
!c  990 CONTINUE
!c      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
!c      IF (PUNCH.EQ.1) CLOSE(7)
!c  999 CONTINUE
return
!xx    3 FORMAT( 1H ,'*** ABIC(',D20.10,'  ) =  ',D20.10 )                 
!xx  600 FORMAT(1H ,'AVABIC =',F10.2)                                      
!xx  602 FORMAT(1H ,12D11.3)                                               
!xx  614 FORMAT( 1H ,'ORIGINAL DATA',/,1H , 20A4 )                         
!xx  615 FORMAT(1H ,'SMOOTHED=TREND+SEASONAL+TRADING.DAY.COMP')            
!xx  630 FORMAT(1H ,'***  LAST ',I3,' VALUES ARE FOCASTED  ***')           
!xx  616 FORMAT( 1H ,'TREND')                                              
!xx  617 FORMAT( 1H ,'SEASONAL')                                           
!xx  618 FORMAT(1H ,'IRREGULAR=ORIGINAL DATA-TREND-SEASONAL-TRADING.',     
!xx     *'DAY.COMP')                                                       
!xx  619 FORMAT(1H ,'ADJUSTED=ORIGINAL DATA-SEASONAL-TRADING.DAY.COMP-OCF')
!xx  622 FORMAT( 1H ,'IRREGULAR ( SCALED BY THE STANDARD DEVIATION )')     
!xx  610 FORMAT(1H ,                                                       
!xx     * 'LOGT  =',I10,/,                                                 
!xx     *' MT    =',I10,/,                                                 
!xx     *' RLIM  =',D10.3,/,                                               
!xx     *' PERIOD=',I10,/,                                                 
!xx     *' SPAN  =',I10,/,                                                 
!xx     *' SHIFT =',I10,/,                                                 
!xx     *' FOCAST=',I10,/,                                                 
!xx     *' ORDER =',I10,/,                                                 
!xx     *' SORDER=',I10,/,                                                 
!xx     *' RIGID =',D10.3,/,                                               
!xx     *' YEAR  =',I10,/,                                                 
!xx     *' MONTH =',I10,/,                                                 
!xx     *' WTRD  =',D10.3,/,                                               
!xx     *' IOUTD =',I10,/,                                                 
!xx     *' SPEC  =',I10,/,                                                 
!xx     *' PUNCH =',I10,/,                                                 
!xx     *' ZERSUM=',D10.3,/,                                               
!xx     *' DELTA =',D10.3)                                                 
!xx  611 FORMAT( 1H  )                                                     
!xx  612 FORMAT(1H ,'D3    =',10F10.6)                                     
!xx  640 FORMAT(1H ,'TRADING-DAY COMPONENT')                               
!xx  641 FORMAT(1H ,'OUTLIER CORRECTION FACTOR')                           
!xx  702 FORMAT(6D12.5)                                                    
!xx  700 FORMAT('(6D12.5)')                                                
!xx  714 FORMAT(20A4)                                                      
!xx  715 FORMAT('SMOOTHED=TREND+SEASONAL+TRADING.DAY.COMP')                
!xx  716 FORMAT('TREND')                                                   
!xx  717 FORMAT('SEASONAL')                                                
!xx  718 FORMAT('IRREGULAR')                                               
!xx  719 FORMAT('ADJUSTED=DATA-SEASONAL-TRADING.DAY.COMP-OCF')             
!xx  720 FORMAT(I5)                                                        
!xx  721 FORMAT('OUTLIER CORRECTION FACTOR')                               
!xx  730 FORMAT('TRADING DAY COMPONENT')                                   
!xx  703 FORMAT('DADJ')                                                    
!xx  777 FORMAT('MISSING OBSERVATION INTERPOLATED DATA')                   
!xx  776 FORMAT(1H ,'MISSING OBSERVATION INTERPOLATED DATA')               
end
subroutine  add(x,mx,y,my,z,mz)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES                                          
!          X = Y + Z.                                                   
!     INPUTS:                                                           
!       X:     MX-VECTOR                                                
!       Y:     MY-VECTOR                                                
!       Z:     MZ-VECTOR                                                
!                                                                       
!xx      IMPLICIT REAL*8 ( A-H,O-Z )                                       
!c      DIMENSION X(1), Y(1), Z(1)                                        
!xx      DIMENSION X(MX), Y(MY), Z(MZ)                                        
integer mx, my, mz
real(dp) x(mx), y(my), z(mz), tem
! local
integer i
do 100 i=1,mx
tem = 0.d0
if( i .le. my )  tem = y(i)
if( i .le. mz )  tem = tem + z(i)
!xx  100 X(I) = TEM
x(i) = tem
100 continue
return
end
subroutine clear(x,m,n,mj,i0,j0)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE CLEARS MATRIX X.                                  
!     INPUTS:                                                           
!       X:     M*N MATRIX                                               
!       MJ:    ABSOLUTE DIMENSION OF X                                  
!       I0:    ABSOLUTE POSITION OF THE FIRST ROW OF X                  
!       J0:    ABSOLUTE POSITION OF THE FIRST COLUMN OF X               
!                                                                       
!xx      IMPLICIT REAL*8 (A-H,O-Z )                                        
!x      DIMENSION X(MJ,1)                                                 
!xx      DIMENSION X(MJ, I0+N-1)                                                 
integer m, n, mj, i0, j0
real(dp) x(mj, i0+n-1)
! local
integer i, i0m1, j, j0m1
!
i0m1 = i0 - 1
j0m1 = j0 - 1
!xx      DO 10 J=1,N
do 20 j=1,n
do 10 i=1,m
!xx   10 X(I0M1+I,J0M1+J) = 0.D0
x(i0m1+i,j0m1+j) = 0.d0
10 continue
20 continue
return
end
!c      SUBROUTINE COPY(X,MX,NX,MMX,IX,JX,Y,MY,NY,MMY,IY,JY)              
!xx      SUBROUTINE BCOPY(X,MX,NX,MMX,IX,JX,Y,MY,NY,MMY,IY,JY)
subroutine bcopy(x,mx,nx,ix,jx,y,my,ny,iy,jy)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COPIES Y INTO X.                                  
!     INPUTS:                                                           
!       X:     MX*NX MATRIX                                             
!       MMX:   ABSOLUTE DIMENSION OF X                                  
!       IX:    ABSOLUTE POSITION OF THE FIRST ROW OF X                  
!       JX:    ABSOLUTE POSITION OF THE FIRST COLUMN OF X               
!       Y:     MY*NY MATRIX                                             
!       MMY:   ABSOLUTE DIMENSION OF Y                                  
!       IY:    ABSOLUTE POSITION OF THE FIRST ROW OF Y                  
!       JY:    ABSOLUTE POSITION OF THE FIRST COLUMN OF Y               
!                                                                       
!xx      IMPLICIT  REAL*8 ( A-H,O-Z )                                      
!x      DIMENSION X(MMX,1), Y(MMY,1)
!xx      DIMENSION X(MMX,JX+NX-1), Y(MMY,JY+NY-1)
!xx      DIMENSION X(MX+IX-1,JX+NX-1), Y(MY+IY-1,NY+JY-1)
integer mx, nx, ix, jx, my, ny, iy, jy
real(dp) x(mx+ix-1,jx+nx-1), y(my+iy-1,ny+jy-1), tem
! local
integer i, ixm1, iym1, j, jxm1, jym1
!
ixm1 = ix-1
jxm1 = jx - 1
iym1 = iy - 1
jym1 = jy - 1
do 100 j=1,nx
!xx      DO 50 I=1,MX
do 60 i=1,mx
tem = 0.d0
if( i .gt. my ) go to 50
if( j .gt. ny ) go to 50
tem = y(iym1+i,jym1+j)
50 x(ixm1+i,jxm1+j) = tem
60 continue
100 continue
return
end
subroutine decode(seas0,trend0,est0,adj0,irreg0,tdc0,w,&
!c     * A,Y,NN,NF,WEEK,ERR,PSDS,PSDT,SQE)                                
&a,y,nn,nf,week,err,psds,psdt,sqe, ip,year,nday)
  use timsac_kinds, only: dp
  implicit none
!   THIS SUBROUTINE COMPUTES                                            
!       TREND0                                                          
!       SEAS0                                                           
!       EST0=TREND0 + SEAS                                              
!       ADJ0=Y - SEAS0                                                  
!       IRREG0=Y - EST0                                                 
!                                                                       
!xx      IMPLICIT REAL*8 (A-H,O-Z )                                        
!xx      INTEGER*4 YEAR                                                    
!xx      REAL*8 IRREG0                                                     
!c      DIMENSION A(1),Y(1),SEAS0(1),TREND0(1),EST0(1),ADJ0(1),IRREG0(1)  
!c     * ,W(1),WEEK(7,1),TDC0(1), PSDT(500),PSDS(500),ERR(1000)           
!x      DIMENSION A(2*(NN+NF)+NDAY+6),Y(1),SEAS0(NN+NF),TREND0(NN+NF),
!x     * EST0(NN+NF),ADJ0(NN+NF),IRREG0(NN+NF),W(8),WEEK(7,1),
!xx      DIMENSION A(2*(NN+NF)+NDAY+6),Y(NN),SEAS0(NN+NF),TREND0(NN+NF),
!xx     * EST0(NN+NF),ADJ0(NN+NF),IRREG0(NN+NF),W(NDAY+6),WEEK(7,NN+NF),
!xx     * TDC0(NN+NF), PSDT(NN+NF),PSDS(NN+NF),ERR(2*(NN+NF))
integer nn, nf, ip, year, nday
real(dp) seas0(nn+nf), trend0(nn+nf), est0(nn+nf),&
&adj0(nn+nf), irreg0(nn+nf), tdc0(nn+nf),&
&w(nday+6), a(2*(nn+nf)+nday+6), y(nn),&
&week(7,nn+nf), err(2*(nn+nf)), psds(nn+nf),&
&psdt(nn+nf), sqe, sd2
!c      COMMON /IDATA/ IP,IDUMMY(2),YEAR,NDAY 
! local
integer i, i1, i2, n, nr, ntem, n7
!                                                                       
n=nn+nf
nr = 2
if( ip .eq. 1 )   nr = 1
call  clear( seas0,n,1,n,1,1 )
call  clear( psds,n,1,n,1,1)
sd2 = dsqrt(sqe) * 2.d0
do 10 i=1,n
i1=nr*(i-1)+1
i2=nr*i
trend0(i)=a(i1)
psdt(i)=dsqrt(err(i1))*sd2
if(ip .le. 1) go to 10
seas0(i)=a(i2)
psds(i)=dsqrt(err(i2))*sd2
10 continue
if(year .eq. 0) go to 20
ntem = n*2+1
n7 = 6 + nday
!c      CALL COPY(W,N7,1,8,1,1,A,N7,1,N,NTEM,1)
!x      CALL PRDCT(TDC0,1,N,1,W,1,N7,1,WEEK,N7,N,N7)     
!x      CALL BCOPY(W,N7,1,8,1,1,A,N7,1,N,NTEM,1)
!xx      CALL BCOPY(W,N7,1,N7,1,1,A,N7,1,N,NTEM,1)
call bcopy(w,n7,1,1,1,a,n7,1,ntem,1)
call prdct(tdc0,1,n,1,w,1,n7,1,week,n7,n,7)
20 continue
call add(est0,n,trend0,n,seas0,n)
if(year .ne. 0) call add(est0,n,est0,n,tdc0,n)
!xx      CALL SBTRCT(ADJ0,N,Y,N,SEAS0,N)
call sbtrct(adj0,n,y,nn,seas0,n)
if(year .ne. 0) call sbtrct(adj0,n,adj0,n,tdc0,n)
!xx      CALL SBTRCT(IRREG0,N,Y,N,EST0,N)                                  
call sbtrct(irreg0,n,y,nn,est0,n)
return
end
!c      SUBROUTINE  HUSHLD( X,N,K,MJ1,ICNT )                              
subroutine  bhushld( x,n,k,mj1,icnt )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!                                                                       
!     THIS SUBROUTINE  TRANSFORMS MATRIX X INTO AN UPPER TRIANGULAR FORM
!     BY HOUSEHOLDER TRANSFORMATION.                                    
!                                                                       
!       INPUTS:                                                         
!          X:     ORIGINAL N*K DATA MATRIX                              
!          MJ1:   ABSOLUTE DIMENSION OF X                               
!          N:  NUMBER OF ROWS OF X, NOT GREATER THAN MJ1                
!          K:     NUMBER OF COLUMNS OF X. NOT GREATER THAN MJ1          
!          ICNT:  =0  WHEN X IS A FULL N*K MATRIX                       
!                 =L (L.NE.0) WHEN X IS COMPOSED OF TWO UPPER TRIANGULAR
!                             MATRICES. L SHOULD BE EQUAL TO THE NUMBER 
!                             OF NON ZERO ROWS OF THE SECOND MATRIX.    
!                             SECOND MATRIX MUST BE ROTATED AND STORED  
!                             AT THE BOTTOM LEFT PART OF THE FIRST      
!                             MATRIX.                                   
!       OUTPUT:                                                         
!          X:     IN UPPER TRIANGULAR FORM                              
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!c      DIMENSION  X(MJ1,1) , D(1000)                                     
!x      DIMENSION  X(MJ1,1) , D(N)                                     
!xx      DIMENSION  X(MJ1,K) , D(N)
integer n, k, mj1, icnt
real(dp) x(mj1,k)
! local
integer i, ii, ii1, ii10, iiotem, iitem, j, jtem, ktem, mnk
real(dp) d(n), tol, diio, h, absld, f, g, s
!                                                                       
tol=1.0d-38
diio=0.0d00
!                                                                       
mnk=k
if(n.le.k) mnk=n-1
!xx      DO 100 II=1,MNK                                                   
do 101 ii=1,mnk
h = 0.0d00
iiotem = ii
iitem = ii
if( icnt .le. 0 )   go to 5
h = x(ii,ii)*x(ii,ii)
iitem = k+1-ii
iiotem = n+1-ii
if( iiotem .le. n-icnt )   iiotem = n-icnt+1
5 continue
do 10  i=iiotem,n
d(i) = x(i,iitem)
absld=dabs(d(i))
if(absld.le.tol) d(i)=0.0d-00
!xx   10       H = H + D(I)*D(I)                                           
h = h + d(i)*d(i)
10 continue
if( h .gt. tol )  go to 20
g = 0.0d00
go to 100
20 g = dsqrt( h )
f=x(ii,ii)
if( f .ge. 0.0d00 )   g = -g
if( icnt .le. 0 )   d(ii) = f-g
if( icnt .gt. 0 )   diio = f-g
h = h - f*g
!                                                                       
!          FORM  (I - D*D'/H) * X, WHERE H = D'D/2                      
!                                                                       
ii1 = iitem+1
ktem = k
if( icnt .le. 0 )   go to 25
ii1 = 1
ktem = iitem-1
25 continue
ii10 = ii1
if( icnt .gt. 0 )   ii10 = iiotem
do 30 i=ii10,n
!xx   30 X(I,IITEM) = 0.D0                                                 
x(i,iitem) = 0.d0
30 continue
if( ii .eq. k )  go to 100
do 60  j=ii1,ktem
s = 0.0d00
jtem = k+1-j
if(icnt .gt. 0 ) s = diio*x(ii,jtem)
do 40  i=iiotem,n
!xx   40 S = S + D(I)*X(I,J)                                               
s = s + d(i)*x(i,j)
40 continue
s = s/h
if(icnt .gt. 0) x(ii,jtem) = x(ii,jtem) - diio*s
do 50  i=iiotem,n
!xx   50 X(I,J) = X(I,J) - D(I)*S
x(i,j) = x(i,j) - d(i)*s
50 continue
60 continue
100 x(ii,ii) = g
101 continue
!                                                                       
return
!                                                                       
end
subroutine  prdct(x,mx,nx,mmx,y,my,ny,mmy,z,mz,nz,mmz)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES                                          
!          X = Y * Z                                                    
!     INPUTS:                                                           
!       X:     MX*NX MATRIX                                             
!       MMX:   ABSOLUTE DIMENSION OF X                                  
!       Y:     MY*NY MATRIX                                             
!       MMY:   ABSOLUTE DIMENSION OF Y                                  
!       Z:     MZ*NZ MATRIX                                             
!       MMZ:   ABSOLUTE DIMENSION OF Z                                  
!                                                                       
!xx      IMPLICIT REAL*8 ( A-H,O-Z )                                       
!x      DIMENSION X(MMX,1), Y(MMY,1), Z(MMZ,1)                            
!xx      DIMENSION X(MMX,NX), Y(MMY,NY), Z(MMZ,NZ)
integer mx, nx, mmx, my, ny, mmy, mz, nz, mmz
real(dp) x(mmx,nx), y(mmy,ny), z(mmz,nz), sum
! local
integer i, j, k, kk
!
kk = ny
if( kk .gt. mz ) kk = mz
do 100 j=1,nx
!xx      DO 50 I=1,MX
do 51 i=1,mx
sum = 0.d0
if( i .gt. my ) go to 50
if( j .gt. nz ) go to 50
do 20 k=1,kk
!xx   20 SUM = SUM + Y(I,K)*Z(K,J)
sum = sum + y(i,k)*z(k,j)
20 continue
50 x(i,j) = sum
51 continue
100 continue
return
end
subroutine  sbtrct(x,mx,y,my,z,mz)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES                                          
!          X = Y - Z                                                    
!     INPUTS:                                                           
!       X:     MX-VECTOR                                                
!       Y:     MY-VECTOR                                                
!       Z:     MZ-VECTOR                                                
!                                                                       
!xx      IMPLICIT REAL*8 ( A-H,O-Z )                                       
!c      DIMENSION X(1), Y(1), Z(1)                                        
!xx      DIMENSION X(MX), Y(MY), Z(MZ)
integer mx, my, mz
real(dp) x(mx), y(my), z(mz), tem
! local
integer i
!
do 100 i=1,mx
tem = 0.d0
if( i .le. my )  tem = y(i)
if( i .le. mz )  tem = tem - z(i)
!xx  100 X(I) = TEM
x(i) = tem
100 continue
return
end
subroutine subsea(abicm,season,trend,est,adj,irreg,tdc,&
!c     *   FSEA,FTRN,YS,N,NF,RLIM,WEEK,DC,IDC,H,NH,F,H2,N2,ITRN,          
!c     *  IARS,ARFS,IART,ARFT,IARN,ARFN,PSDT,PSDS)                        
&fsea,lfsea,ftrn,ys,n,nf,rlim,week,idc,nh,f,n2,itrn,iars,arfs,&
&iart,arft,iarn,arfn,psdt,psds,npf,period,iord,isod,year,nday,&
&logt,alpha,beta,gamma,zer,smth,smth2,dd,wtrd,delta)
  use timsac_kinds, only: dp
  implicit none
!     SEASONAL DECOMPOSITION PROCEDURE                                  
!     FOR THE DEFINITIONS OF THE VARIABLES APPEARING IN THE ARGUMENTS,  
!     SEE THE COMMENTS IN THE MAIN ROUTINE                              
!xx      IMPLICIT  REAL*8  ( A-H,O-Z )                                    
!c      REAL*8  IRREG0, IRREG                                             
!xx      REAL*8  IRREG                                             
!xx      INTEGER*4 PERIOD,YEAR                                             
!c      DIMENSION   FSEA(1), FTRN(1), DC(IDC,1), YS(1), H2(N2,1)
!c     *   ,ERR(1000),PSDT(500),PSDS(500)                                 
!c     *  ,SEASON(1),TREND(1),EST(1),ADJ(1),IRREG(1),A(1000)
!c     * ,WEEK(7,1),TDC(1),WEEK0(8),WEEK1(8),H(NH,1),F(1)                 
!c      DIMENSION  DTRN(500),DSEAS(500), ARFS(1), ARFT(1), ARFN(1)        
!xx      DIMENSION  SEASON(NPF),TREND(NPF),EST(NPF),ADJ(NPF),
!xx     *   IRREG(NPF),TDC(NPF),
!x     *   FSEA(LFSEA),FTRN(IORD+3),YS(N),WEEK(7,1),F(NPF+1),ARFS(1),
!x     *   ARFT(3),ARFN(1),PSDT(NPF),PSDS(NPF)
!xx     *   FSEA(LFSEA),FTRN(IORD+3),YS(N),WEEK(7,1),F(NPF+1),ARFS(3),
!xx     *   ARFT(3),ARFN(3),PSDT(NPF),PSDS(NPF)
!x      DIMENSION  DC(IDC,2*NPF+N2),H(NH,NPF),H2(N2,2*NPF+N2),WEEK0(8),
!x     *   WEEK1(8),ERR(2*(N+NF)+NDAY+7),A(2*(N+NF)+NDAY+7),DTRN(NPF),
!xx      DIMENSION  DC(IDC,2*NPF+N2),H(NH,NPF),H2(N2,2*NPF+N2),WEEK0(7),
!xx     *   WEEK1(7),ERR(2*(N+NF)+NDAY+7),A(2*(N+NF)+NDAY+7),DTRN(NPF),
!xx     *  DSEAS(NPF)
integer lfsea, n, nf, idc, nh, n2, itrn, iars, iart, iarn, npf,&
&period, iord, isod, year, nday, logt
real(dp) abicm, season(npf), trend(npf), est(npf),&
&adj(npf), irreg(npf), tdc(npf), fsea(lfsea),&
&ftrn(iord+3), ys(n), rlim, week(7,1),  f(npf+1),&
&arfs(3), arft(3), arfn(3), psdt(npf), psds(npf),&
&alpha, beta, gamma, zer, smth, smth2, dd, wtrd,&
&delta
! local
integer i, iiii, icount, iflag, itrn0, ipm1, j, k, m1, mode, nd,&
&n7, n2m1, ndtem, nmj
real(dp) dc(idc,2*npf+n2), h(nh,npf),h2(n2,2*npf+n2),&
&week0(7), week1(7), err(2*(n+nf)+nday+7),&
&a(2*(n+nf)+nday+7), dtrn(npf), dseas(npf), dmax0,&
&dmin0, ro, dmin, an, ann, alndtd, alndt0, wt,&
&tem, alndn, alsqe, sqe, abic, ajacob, aprd, sstr,&
&ssea, ssir, ssas, sas
!c      COMMON /IDATA/ PERIOD,IORD,ISOD,YEAR,NDAY,IFIX                    
!c      COMMON /RDATA/ ALPHA,BETA,GAMMA,ZER,SMTH,SMTH2,DD,WTRD,DELTA      
!c      COMMON /ILOGT/ LOGT,ISHRNK                                        
iflag=0
dmax0 = 1000.d0
dmin0 = 1.0d0
mode = 0
ro = 1.41421d0
if(itrn .ne. 0) ro = dsqrt(ro)
nd=(n+nf)*2
if(period .eq. 1) nd=n+nf
n7 = nday + 6
if(year .ne. 0) nd=nd+n7
abicm = 1.d50
!-----
dmin = dmin0
ann = nd
!-----                                                     
!                                                                       
alndt0=0.d0
itrn0 = 30
dd=dmin
if(itrn.eq.1) dd=5.d0
do 9999  iiii=1,itrn0
!                                                                       
!     BASIC ROUTINE : SEASONAL ADJUSTMENT UNDER GIVEN PRIOR DISTRIBUTION
!                                                                       
!                                                                       
!      CONSTRUCTION AND HOUSEHOLDER TRANSFORMATION OF MATRIX DC         
!                                                                       
wt=dd*smth
if(iiii.gt.1)go to 8888
if(period.eq.1)go to 8888
!c      CALL SETDC(H,NH,F,M1,FSEA,N+NF,SMTH2,ZER,IARS,ARFS,IARN,ARFN)     
!c     *        ALPHA,BETA,GAMMA,WTRD,DELTA,PERIOD,IORD,ISOD,YEAR)
call setdc(h,nh,f,m1,fsea,n+nf,smth2,zer,iars,arfs,iarn,arfn,&
!x     *        ALPHA,BETA,GAMMA,WTRD,DELTA,PERIOD,ISOD,YEAR)
!xx     *        ALPHA,BETA,GAMMA,WTRD,DELTA,PERIOD,ISOD,YEAR,NPF)
&beta,gamma,period,isod,npf)
alndt0=0.d0
do 2233 i=1,m1
tem=dabs(h(1,i))
!xx 2233 ALNDT0=ALNDT0+DLOG(TEM)                                           
alndt0=alndt0+dlog(tem)
2233 continue
8888 continue
!                                                                       
!     CALCULATION OF LOG(DET(DC'DC))*0.5                                
!                                                                       
!     --------------------                                              

alndtd=(n+nf)*dlog(wt)+iord*dlog(alpha)+alndt0
if(period.ne.1)alndtd=alndtd+(n+nf)*dlog(dd)
!     --------------------                                              
!                                                                       
!     CONSTRUCTION AND HOUSEHOLDER TRANSFORMATION OF MATRIX DCX         
!                                                                       
!xx 4567 CONTINUE                                                          
!c      CALL SETX(DC,IDC,H2,N2,M1,ICOUNT,FTRN,N+NF,H,NH,WT                
!c     *             ,YS,N,RLIM,WEEK,N7,ALNDTD,F,DD,IART,ARFT)            
call setx(dc,idc,h2,n2,m1,icount,ftrn,n+nf,h,nh,wt,&
&ys,n,rlim,week,n7,alndtd,f,dd,iart,arft,alpha,&
!xx     *  BETA,GAMMA,WTRD,DELTA,PERIOD,IORD,ISOD,YEAR,NPF)
&wtrd,delta,period,iord,year,npf)
!                                                                       
!     LEAST SQUARES COMPUTATION                                         
!                                                                       
k=m1 + n2
sqe = h2(n2,k)**2
!                                                                       
!     INTERPRETATION OF THE SOLUTION                                    
!                                                                       
!     ****                                                              
!                                                                       
!                                                                       
!     ABIC COMPUTATION                                                  
!                                                                       
an = icount
ann = icount + nd
alndn=0.d0
do 3344 i=1,m1
tem=dabs(dc(1,i))
!xx 3344 ALNDN=ALNDN + DLOG(TEM)                                           
alndn=alndn + dlog(tem)
3344 continue
if(n2 .eq. 1) go to 3346
n2m1=n2-1
do 3345 i=1,n2m1
tem=dabs(h2(i,m1+i))
!xx 3345 ALNDN=ALNDN+DLOG(TEM)                                             
alndn=alndn+dlog(tem)
3345 continue
3346 continue
alsqe=an*dlog(sqe/an)
abic=alsqe + 2.d0*(alndn-alndtd)
if(year .ne. 0 .and. wtrd .le. 0.d0)abic=abic+n7*2.d0
!c      WRITE( 6,3 )    DD, ABIC, ALSQE, ALNDN, ALNDTD                    
!                                                                       
!     END OF BASIC ROUTINE                                              
!                                                                       
!                                                                       
!     MINIMUM ABIC PROCEDURE                                            
!                                                                       
if(abic.ge.abicm) go to 9000
if(abicm-abic .lt. 0.0001d0) go to 9000
abicm = abic
dmin = dd
dd=ro*dd
if(iiii .eq. 2) mode = 1
!xx 5678 CONTINUE                                                          
go to 2345
9000 if(mode .eq. 1) go to 6000
9001 mode = 1
ro = 1.d0/ro
dd = dd*ro*ro
2345 if(dd .le. dmax0) go to 1234
if(mode .eq. 0) go to 9001
dd = dmax0
iflag=1
go to 6000
1234 continue
if(dd .ge. dmin0) go to 9999
dd = dmin0
iflag=-1
go to 6000
9999 continue
!                                                                       
6000 continue
dd=dmin
wt = dd*smth
!c      CALL SETX(DC,IDC,H2,N2,M1,ICOUNT,FTRN,N+NF,H,NH,WT                
!c     *             ,YS,N,RLIM,WEEK,N7,ALNDTD,F,DD,IART,ARFT)            
call setx(dc,idc,h2,n2,m1,icount,ftrn,n+nf,h,nh,wt,&
&ys,n,rlim,week,n7,alndtd,f,dd,iart,arft,alpha,&
!xx     *   BETA,GAMMA,WTRD,DELTA,PERIOD,IORD,ISOD,YEAR,NPF)
&wtrd,delta,period,iord,year,npf)
ndtem=nd + 1
!c      CALL SOLVE(DC,IDC,H2,N2,A,M1,SQE,NDTEM,ERR)                       
call bsolve(dc,idc,h2,n2,a,m1,sqe,ndtem,err)
sqe=sqe/ann
call decode(season,trend,est,adj,irreg,tdc,week0,&
!c     *            A,YS,N,NF,WEEK,ERR,PSDS,PSDT,SQE)                     
&a,ys,n,nf,week,err,psds,psdt,sqe,period,year,nday)
if(logt .eq. 0) go to 6200
ajacob=0.d0
do 6100 i=1,n
!xx 6100 IF(YS(I) .LT. RLIM .OR. RLIM .LE. 0.D0) AJACOB=AJACOB+YS(I)       
if(ys(i) .lt. rlim .or. rlim .le. 0.d0) ajacob=ajacob+ys(i)
6100 continue
ajacob=ajacob+ajacob
abicm=abicm+ajacob
6200 continue
do 6300 i=1,7
week1(i) = week0(i)
if(logt .ne. 0 .and. year .ne. 0) week1(i) = dexp(week1(i))
6300 continue
!xx    3 FORMAT( 1H ,'    ABIC(',D20.10,'  ) =  ',F10.2,                   
!c     *3X,'ALSQE=',D13.5,3X,'ALNDN=',D13.5,3X,'ALNDTD='D13.5)            
!xx     *3X,'ALSQE=',D13.5,3X,'ALNDN=',D13.5,3X,'ALNDTD=',D13.5)           
do 3320 i=1,n
dtrn(i)=trend(i)
!xx 3320 DSEAS(I)=SEASON(I)
dseas(i)=season(i)
3320 continue
if(iord .eq. 0) go to 3323
!xx      DO 3321 J=1,IORD
do 3322 j=1,iord
nmj=n-j
do 3321 i=1,nmj
!xx 3321 DTRN(I)=DTRN(I+1)-DTRN(I)
dtrn(i)=dtrn(i+1)-dtrn(i)
3321 continue
3322 continue
3323 if(isod .eq. 0) go to 3325
!xx      DO 3324 J=1,ISOD                                                  
do 3327 j=1,isod
nmj=n-j*period
do 3324 i=1,nmj
!xx 3324 DSEAS(I)=DSEAS(I+PERIOD)-DSEAS(I)
dseas(i)=dseas(i+period)-dseas(i)
3324 continue
3327 continue
3325 sstr=0.d0
ssea=0.d0
ssir=0.d0
ssas=0.d0
sas=0.d0
ipm1=period-1
do 3326 i=1,ipm1
!xx 3326 SAS=SAS+SEASON(I)
sas=sas+season(i)
3326 continue
do 3400 i=1,n
if(rlim .gt. 0.d0 .and. irreg(i) .gt. rlim) go to 3400
ssir=ssir+irreg(i)**2
3400 continue
nmj=n-iord
do 3410 i=1,nmj
!xx 3410 SSTR=SSTR+DTRN(I)**2                                              
sstr=sstr+dtrn(i)**2
3410 continue
nmj=n-isod*period
do 3420 i=1,nmj
sas=sas+season(i+ipm1)
ssas=ssas+sas**2
sas=sas-season(i)
!xx 3420 SSEA=SSEA+DSEAS(I)**2                                             
ssea=ssea+dseas(i)**2
3420 continue
aprd=period
ssas=ssas/aprd
!c      WRITE(6,3450) SSIR,SSTR,SSEA,SSAS                                 
!xx 3450 FORMAT(1H ,'SS IRREGULAR =',D13.5,5X,'SS TREND =',D13.5,5X,'SS    
!xx     * SEASONAL =',D13.5,5X,'SS AVSEAS =',D13.5,/,1H )                  
!c      IF(YEAR .NE. 0) WRITE(6,602)  (WEEK1(I),I=1,7)                    
!c      WRITE(6,606)  ABICM, DMIN                                         
!xx  606 FORMAT(1H ,'MINIMUM ABIC =',F10.2,'  ATTAINED AT D =',D20.10 )    
!c      IF(IFLAG .EQ. 1) WRITE(6,600) DMAX0                               
!c      IF(IFLAG .EQ. -1) WRITE(6,601) DMIN0                              
!xx 9876 RETURN                                                            
return
!xx  600 FORMAT(1H ,'**** D IS HITTING THE UPPER BOUND ',D13.5,' ----- TRY'
!xx     * ,  ' LOWER VALUES OF ORDER AND/OR SORDER')                       
!xx  601 FORMAT(1H ,'**** D IS HITTING THE LOWER BOUND ',D13.5,' ----- TRY'
!xx     * ,   ' HIGHER VALUES OF ORDER AND/OR SORDER')                     
!xx  602 FORMAT(1H ,4X,'MON',9X,'TUE',9X,'WED',9X,'THU',9X,                
!xx     *  'FRI',9X,'SAT',9X,'SUN',/,1H ,7D12.4)                           
!xx  603 FORMAT(1H ,'SHRINKAGE FACTORS ARE  ',D12.5,' FOR TREND, ',        
!xx     *  D12.5,' FOR SEASONAL.')                                         
end
subroutine calend(week,year0,month0,n)
  use timsac_kinds, only: dp
  implicit none
!      THIS SUBROUTINE COMPUTES THE DAYS-OF-WEEK DISTRIBUTION OF        
!     N SUCCESSIVE MONTHS STARTING AT MONTH0 OF YEAR0                   
!     NOTE:  THIS SUBROUTINE WORKS FOR YEARS                            
!              AD.1901 - AD.2099                                        
!                                                                       
!xx      IMPLICIT INTEGER*4 (A-Z)                                          
!x      REAL*8 WEEK(7,1),W0(8)                                            
!xx      REAL*8 WEEK(7,N),W0(8)
integer year0, month0, n
real(dp) week(7,n)
! local
integer i, j, dyear, year, leap, y, l, d, day, month, diff, wday
real(dp) w0(8)
!                                                                       
dyear=(month0-1)/12
if(month0 .ge. 1) go to 20
dyear=-month0
dyear=dyear/12+1
dyear=-dyear
20 continue
month=month0-12*dyear
year=year0+dyear
leap=mod(year,4)
!  DAY-OF-WEEK OF THE FIRST DAY OF 'YEAR'                               
y = year-1901
l=y/4
d=y+l+2
day=mod(d,7) + 1
!  DAY-OF-WEEK OF THE FIRST DAY OF 'MONTH' OF 'YEAR'                    
!xx      GO TO (200,203,203,206,201,204,206,202,205,200,203,205), MONTH    
!xx  201 DAY=DAY+1                                                         
!xx      GO TO 200                                                         
!xx  202 DAY=DAY+2                                                         
!xx      GO TO 200                                                         
!xx  203 DAY=DAY+3                                                         
!xx      GO TO 200                                                         
!xx  204 DAY=DAY+4                                                         
!xx      GO TO 200                                                         
!xx  205 DAY=DAY+5                                                         
!xx      GO TO 200                                                         
!xx  206 DAY=DAY+6                                                         
if (month .eq. 1 .or. month .eq. 10) go to 200
if (month .eq. 5) then
day=day+1
else if (month .eq. 8) then
day=day+2
else if (month .eq. 6) then
day=day+4
else if (month .eq. 9 .or. month .eq. 12) then
day=day+5
else if (month .eq. 7 .or. month .eq. 4) then
day=day+6
else
!       IF(MONTH .EQ. 2 .OR. MONTH .EQ. 3 .OR. MONTH .EQ. 11)
day=day+3
end if
200 if(leap .eq. 0 .and. month .ge. 3) day=day+1
if(day .gt. 7) day=day-7
!  ITERATION                                                            
do 100 i=1,n
!xx      DO 10 J=1,7                                                       
!xx   10 W0(J) = 4.D0
w0(1:7) = 4.d0
!xx      GO TO (331,328,331,330,331,330,331,331,330,331,330,331), MONTH    
if (month .eq. 2) go to 328
if (month .eq. 4 .or. month .eq. 6 .or.&
&month .eq. 9 .or. month .eq. 11) go to 330
!    IF (MONTH .EQ. 1 .OR. MONTH .EQ. 3 .OR. MONTH .EQ. 5 .OR. MONTH .EQ. 7
!        .OR. MONTH .EQ. 8 .OR. MONTH .EQ. 10 .OR. MONTH .EQ. 12) GO TO 331
!xx  331 DIFF=3                                                            
diff=3
w0(8)=31.d0
go to 300
330 diff=2
w0(8)=30.d0
go to 300
328 diff=0
w0(8)=28.d0
if(leap .ne. 0) go to 50
diff=1
w0(8)=29.d0
300 wday=8-day
do 400 j=1,diff
w0(wday) = 5.d0
if(j .eq. diff) go to 50
wday=wday-1
!xx  400 IF(WDAY .EQ. 0) WDAY=7
if(wday .eq. 0) wday=7
400 continue
50 continue
do 410 j=1,7
!xx  410 WEEK(J,I)=W0(J)-30.4375D0/7.D0
week(j,i)=w0(j)-30.4375d0/7.d0
410 continue
if(i .eq. n) go to 900
day = day + diff
if(day .gt. 7) day = day - 7
month=month+1
if(month .le. 12) go to 100
month=1
year=year+1
leap=mod(year,4)
100 continue
900 return
end
!                                                                       
!C      FUNCTION AMAX(A,N)                                                
real(dp) function amax(a,n)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE                                                 
!     MAXIMUM OF A(I)(I=1,N) SEARCH                                     
integer n
real(dp) a(n)
! local
integer i
!
amax=a(1)
do 10 i=2,n
if(amax.lt.a(i)) amax=a(i)
10 continue
return
end
!                                                                       
!C      FUNCTION AMIN(A,N)                                                
real(dp) function amin(a,n)
  use timsac_kinds, only: dp
  implicit none
!     COMMON SUBROUTINE                                                 
!     MINIMUM OF A(I)(I=1,N) SEARCH                                     
integer n
real(dp) a(n)
! local
integer i
!
amin=a(1)
do 10 i=2,n
if(amin.gt.a(i)) amin=a(i)
10 continue
return
end
!                                                                       
subroutine  poolav( z,k,x,sd1 )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE SEARCHES FOR THE MINIMUM OF                       
!          F(X;Z) = (Z(1)-X(1))**2 + ... + (Z(K)-X(K))**2               
!     SUBJECT TO X(1)<X(2)< ... < X(K), BY THE POOL-ADJACENT-VIOLALORS  
!     ALGORITHM.                                                        
!                                                                       
!     INPUTS:                                                           
!        (Z(I),I=1,K): DATA                                             
!        K:            NUMBER OF DATA                                   
!     OUTPUTS:                                                          
!        (X(I),I=1,K): VECTOR OF MINIMIZING SOLUTION                    
!        SD1:          MINIMUM OF F(X;Z)                                
!                                                                       
!xx      IMPLICIT  REAL*8( A-H,O-Z )                                       
!c      DIMENSION  Z(K), X(K), Y(20)                                      
!xx      DIMENSION  Z(K), X(K), Y(K)                                      
integer k
real(dp) z(k), x(k), sd1
! local
integer i, i0, ifg, j, n0
real(dp) y(k), sum
do 10  i=1,k
!xx   10 X(I) = Z(I)                                                       
x(i) = z(i)
10 continue
!                                                                       
100 do 20  i=2,k
!xx   20 IF( X(I-1) .GT. X(I) )  GO TO 30
if( x(i-1) .gt. x(i) )  go to 30
20 continue
go to 300
!                                                                       
30 ifg = 0
do 40  i=1,k
!xx   40 Y(I) = X(I)                                                       
y(i) = x(i)
40 continue
!
n0=1
do 200  i=1,k-1
i0 = i
!                                                                       
if( x(i) .lt. x(i+1) )  go to 110
if( i .eq. k-1  .and.  ifg .ne. 0 )  go to 90
if( ifg .ne. 0 )  go to 200
ifg = 1
n0 = i
if( i .eq. k-1 )  go to 90
go to 200
!C   90 I = K
90 i0 = k
go to 115
!                                                                       
110 if( ifg .eq. 0 )  go to 200
ifg = 0
115 continue
sum = 0.d0
!C      DO 120  J=N0,I
do 120 j=n0,i0
!xx  120 SUM = SUM + Y(J)                                                  
sum = sum + y(j)
120 continue
!C      SUM = SUM/(I-N0+1)
sum = sum/(i0-n0+1)
!C      DO 130  J=N0,I
do 130  j=n0,i0
!xx  130 Y(J) = SUM                                                        
y(j) = sum
130 continue
!                                                                       
200 continue
!                                                                       
do 210  i=1,k
!xx  210 X(I) = Y(I)                                                       
x(i) = y(i)
210 continue
go to 100
300 sd1 = 0.d0
do 310  i=1,k
!xx  310 SD1 = SD1 + (X(I) - Z(I))**2                                      
sd1 = sd1 + (x(i) - z(i))**2
310 continue
!xx    2 FORMAT( 1H ,10F13.5 )                                             
return
end
subroutine  lkout1( x,n,ind,jsw,f,w )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE COMPUTES THE LIKELIHOOD OF THE MODEL THAT (X(I);  
!     IND(I)=1) ARE THE OUTLIERS.  (MEAN SLIPPAGE TYPE MODEL)           
!                                                                       
!     INPUTS:                                                           
!        (X(I),I=1,N): OBSERVATIONS                                     
!        N:            NUMBER OF OBSERVATIONS                           
!        (IND(I),I=1,N): = 0 ; IF X(I) IS A "NORMAL" OBSERVATION.       
!                        = 1 ; IF X(I) IS AN OUTLIER.                   
!        JSW:   =0;    WHOLE MODELS ARE EVALUATED                       
!               =1;    ONLY THE NATURALLY ORDERED MODEL IS EVALUATED    
!                      (SIMPLIFIED ALGORITHM)                           
!                                                                       
!     OUTPUTS:                                                          
!        F:            LIKELIHOOD OF THE MODEL                          
!        W:                                                             
!                                                                       
!xx      IMPLICIT REAL*8  ( A-H,O-Z )                                      
!c      DIMENSION  X(N), Y(10), Z(10), ZE(10), IND(N), JND(10)            
!xx      DIMENSION  X(N), Y(N), Z(N), ZE(N), IND(N), JND(N)            
integer n, ind(n), jsw
real(dp) x(n), f, w
! local
integer i, ifg, j, k, l, jnd(n)
real(dp) y(n), z(n), ze(n), sum, xmean, sig2, sd
!                                                                       
l = 0
sum = 0.d0
do 10  i=1,n
if( ind(i) .eq. 1 )   go to 10
l = l + 1
sum = sum+x(i)
10 continue
!c      XMEAN = SUM/DFLOAT(L)                                             
xmean = sum/dble(l)
k = n-l
!                                                                       
sum = 0.d0
do 20  i=1,n
!xx   20 IF( IND(I) .EQ. 0 )   SUM = SUM+(X(I)-XMEAN)**2                   
if( ind(i) .eq. 0 )   sum = sum+(x(i)-xmean)**2
20 continue
sig2 = sum/n
w = 1.d0
f = -0.5d0*n*dlog(sig2)
!                                                                       
if( jsw .eq.1 )   return
if( k .le. 1 )   return
!                                                                       
j = 0
do 30  i=1,n
if( ind(i) .eq. 0 )   go to 30
j = j+1
y(j) = x(i)
30 continue
!                                                                       
w = 0.d0
do 40  i=1,k
!xx   40 JND(I) = I                                                        
jnd(i) = i
40 continue
!                                                                       
50 do 60  i=1,k
j = jnd(i)
!xx   60 Z(I) = Y(J)                                                       
z(i) = y(j)
60 continue
call  poolav( z,k,ze,sd )
!                                                                       
w = w + 1.d0/dsqrt(1.d0+sd/sum)**n
!                                                                       
call  permut( jnd,k,ifg )
if( ifg .eq. 0 )   go to 50
!                                                                       
return
!xx  600 FORMAT( 1H ,'IND',40I3 )                                          
!xx  601 FORMAT( 1H ,'F =',D13.5,5X,'FSUM =',D13.5 )                       
end
!c      SUBROUTINE  PRPOST( POST,X,IND,JND,KND,IC,N,L)                    
subroutine  prpost( post,x,ind,jnd,knd,ic,n,l)
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE ARRANGES POST(I), JND(I) AND KND(I) (I=1,IC) IN   
!     DECREASING ORDER OF POST(I), AND DRAWS POSTERIOR PROBABILITY AND  
!     THE ASSUMED OUTLIERS OF THE MODEL WITH POSTERIOR PROBABILITY      
!     GREATER THAN EPS.                                                 
!                                                                       
!     INPUTS:                                                           
!        (POST(I),I=1,IC):   POSTERIOR PROBABILITIES OF THE MODELS      
!        (X(I),I=1,N):       ORIGINAL DATA                              
!        (IND(I),I=1,N):     WORK AREA                                  
!        (JND(I),I=1,IC):    SPECIFICATION OF THE OUTLIERS IN LOW SIDE  
!                            (CODED IN DECIMAL)                         
!        (KND(I),I=1,IC):    SPECIFICATION OF THE OUTLIERS IN HIGH SIDE 
!                            (CODED IN DECIMAL)                         
!        IC:                 NUMBER OF RECORDED MODELS                  
!        N:                  NUMBER OF ORIGINAL DATA                    
!        L:                  NUMBER OF POSSIBLE OUTLIERS IN BOTH SIDES  
!        EPS:                LOWEST LIMIT OF POSTERIOR PROBABILITY TO BE
!                            PRINTED                                    
!                                                                       
!xx      REAL*8  POST, X(N)                                                
!c      DIMENSION  POST(IC), JND(IC),KND(IC), IND(N), Y(10)               
!xx      DIMENSION  POST(IC), JND(IC),KND(IC), IND(N), Y(N)
integer ic, n, l, ind(n), jnd(ic), knd(ic)
real(dp) post(ic), x(n)
! local
integer i, ic1, id, imax, j, jj, kk, nml1
real(dp) y(n), pmax
!c      COMMON /CSPRSS/ ISPRSS                                            
!                                                                       
do 20  i=1,ic
imax = i
pmax = post(i)
do 10  j=i,ic
if( post(j) .le. pmax )   go to 10
imax = j
pmax = post(j)
10 continue
if( imax .eq. i )   go to 20
post(imax) = post(i)
post(i) = pmax
jj = jnd(i)
kk = knd(i)
jnd(i) = jnd(imax)
knd(i) = knd(imax)
jnd(imax) = jj
knd(imax) = kk
20 continue
!xx   30 IC1 = IC                                                          
ic1 = ic
nml1 = n-l+1
do 40  i=1,n
!xx   40 IND(I) = 0                                                        
ind(i) = 0
40 continue
!                                                                       
!c      IF(ISPRSS .EQ. 0) WRITE( 6,4 )                                    
do 100  j=1,ic1
call  binary( jnd(j),l,ind )
call  binary( knd(j),l,ind(nml1) )
id = 0
do 50  i=1,n
if( ind(i) .eq. 0 )    go to 50
id = id+1
y(id) = x(i)
50 continue
!c      IF(ISPRSS .NE. 0) GO TO 100                                       
!c      IF( ID .GE. 1 )   WRITE( 6,5 )   J, POST(J), (Y(I),I=1,ID)        
!c      IF( ID .EQ. 0 )   WRITE( 6,6 )   J, POST(J)                       
100 continue
return
!xx    4 FORMAT( 1H ,10X,'POSTERIOR',10X,'OUTLIERS' )                      
!xx    5 FORMAT( 1H ,I5,F13.5,5X,10F10.3 )                                 
!xx    6 FORMAT( 1H ,I5,F13.5,9X,'NONE' )                                  
end
subroutine permut( ind,k,ifg )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE SEQUENTIALLY SPECIFIES K] CONFIGURATIONS (IND(1), 
!     ...,IND(K)) OBTAINED BY PERMUTING (1,...,K)                       
!                                                                       
!     INPUTS:                                                           
!        (IND(I),I=1,K): FORMER CONFIGURATION                           
!        K:              NUMBER OF ELEMENTS TO BE PERMUTED              
!                                                                       
!     OUTPUTS:                                                          
!        (IND(I),I=1,K): NEW CONFIGURATION                              
!        IFG:            = 0 ; IF THE NEW CONFIGURATION IS OBTAINED     
!                        = 1 ; SEARCH FOR THE CONFIGURATION COMPLETED.  
!                                                                       
!xx      DIMENSION  IND(K)
integer  k, ifg, ind(k)
! local
integer  i, i0, i1, i2, i2m1, imax
!                                                                       
i1 = 1
i2 = 2
i0 = 1
imax = ind(1)
ifg = 0
!                                                                       
10 if( ind(i1) .lt. ind(i2) )   go to 100
if( i2 .eq. i1+1 )   go to 20
i1 = i1+1
go to 10
!                                                                       
20 i2 = i2+1
if( i2 .gt. k )   go to 200
!                                                                       
i2m1 = i2-1
do 30  i=1,i2m1
!xx   30 IF( IND(I) .LE. IND(I2) )   GO TO 40                              
if( ind(i) .le. ind(i2) )   go to 40
30 continue
go to 20
!                                                                       
40 imax = 0
do 50  i=1,i2m1
if( ind(i) .gt. ind(i2) )   go to 50
if( ind(i) .lt. imax )   go to 50
imax = ind(i)
i0 = i
50 continue
!                                                                       
100 ind(i0) = ind(i2)
ind(i2) = imax
i2m1 = i2-1
if( i2 .gt. 2 )   call  isort( ind,i2m1 )
return
!                                                                       
200 ifg = 1
return
!                                                                       
end
subroutine  isort( ind,n )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!     THIS SUBROUTINE ARRANGES IND(I) (I=1,N) IN ORDER OF INCREASING    
!     MAGNITUDE                                                         
!                                                                       
!     INPUTS:                                                           
!        (IND(I),I=1,N): ORIGINAL DATA                                  
!        N:              NUMBER OF DATA                                 
!                                                                       
!     OUTPUT:                                                           
!        (IND(I),I=1,N): REORDERED DATA                                 
!                                                                       
!xx      DIMENSION  IND(N)                                                 
integer n, ind(n)
! local
integer i, ii, imin, j, mini, nm1
!                                                                       
nm1 = n-1
do 20  ii=1,nm1
mini = ind(ii)
imin = ii
do 10  i=ii,n
if( mini .le. ind(i) )   go to 10
mini = ind(i)
imin = i
10 continue
if( imin .eq. ii )   go to 20
j = ind(ii)
ind(ii) = mini
ind(imin) = j
20 continue
return
end
subroutine  binary( m,k,mb )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!       DECIMAL TO BINARY CONVERSION                                    
!                                                                       
!       INPUTS:                                                         
!          M:     NUMBER IN DECIMAL REPRESENTATION                      
!          K:     NUMBER OF BITS USED FOR THE BINARY REPRESENTATION     
!                                                                       
!       OUTPUT:                                                         
!          MB:    NUMBER IN BINARY REPRESENTATION                       
!                                                                       
!x      DIMENSION  MB(1)
!xx      DIMENSION  MB(K)
integer m, k, mb(k)
! local
integer i, l, n
!                                                                       
n = m
do 10  i=1,k
l = n / 2
mb(i) = n - l*2
!xx   10 N = L
n = l
10 continue
return
!                                                                       
end
!c      SUBROUTINE  SRTMIN( X,N,IX )                                      
subroutine  bsrtmin( x,n,ix )
  use timsac_kinds, only: dp
  implicit none
!                                                                       
!       THIS SUBROUTINE ARRANGES X(I) (I=1,N) IN ORDER OF INCREASING    
!       MAGNITUDE OF X(I)                                               
!                                                                       
!       INPUTS:                                                         
!          X:   VECTOR                                                  
!          N:   DIMENSION OF THE VECTOR                                 
!       OUTPUTS:                                                        
!          X:   ARRANGED VECTOR                                         
!          IND: INDEX OF ARRANGED VECTOR                                
!                                                                       
!xx      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
!x      DIMENSION  X(1) , IX(1)                                           
!xx      DIMENSION  X(N) , IX(N)
integer n, ix(n)
real(dp) x(n), xmin, xt
! local
integer i, ii, it, min, nm1
!                                                                       
nm1 = n - 1
do  20     ii=1,nm1
xmin = x(ii)
min = ii
do  10     i=ii,n
if( xmin .lt. x(i) )     go to 10
xmin = x(i)
min = i
10 continue
if( xmin .eq. x(ii) )     go to 20
xt = x(ii)
x(ii) = x(min)
x(min) = xt
it = ix(ii)
ix(ii) = ix(min)
ix(min) = it
20 continue
!                                                                       
return
end
!c      SUBROUTINE OUTLIR( Z,NN,K,ISW,JSW,Y,RLIM )                        
subroutine outlir( z,nn,k,isw,jsw,y,rlim,ioutd,rout )
  use timsac_kinds, only: dp
  implicit none
!     REVISED  MAY 30, 1980                                             
!                                                                       
!     INPUTS:                                                           
!        (Z(I),I=1,N):   ORIGINAL DATA                                  
!        NN:             NUMBER OF DATA                                 
!        K:              MAXIMUM OF THE NUMBER OF OUTLIERS              
!        L:              NUMBER OF POSSIBLE OUTLIERS IN BOTH SIDES.     
!        ISW   =0:       MODIFIED DATA ARE NOT GIVEN                    
!              =2:       THE OBSERVATIONS JUDGED AS OUTLIERS ARE        
!                        REPLACED BY A CONSTANT.                        
!        JSW   =0:       WHOLE MODELS ARE EVALUATED                     
!              =1:       ONLY THE NATURALLY ORDERED MODEL IS EVALUATED  
!                        (SIMPLIFIED ALGORITHM)                         
!        RLIM:           ORIGINAL DATA WHOSE VALUES ARE GREATER         
!                        THAN OR EQUAL TO RLIM ARE TREATEDAS            
!                        MISSING OBSERVATIONS.                          
!     OUTPUTS:                                                          
!        (IX(I),I=1,N):  SUBSCRIPT INDICATING THE ORDER OF THE MAGNITUDE
!                        OF ORIGINAL DATA.                              
!                        I.E.,  Z(IX(1))<Z(IX(2))< ... <Z(IX(N)).       
!                        NOTE THAT Z(IX(I))=X(I).                       
!        (PM(I),I=1,N):  MARGINAL POSTRIOR PROBABILITY THAT X(I) IS AN O
!        (POST(I),I=1,IC):   POSTERIOR PROBABILITIES OF THE MODELS      
!        (JND(I),I=1,IC):    SPECIFICATION OF THE OUTLIERS IN LOW SIDE  
!                            (CODED IN DECIMAL)                         
!        (KND(I),I=1,IC):    SPECIFICATION OF THE OUTLIERS IN HIGH SIDE 
!                            (CODED IN DECIMAL)                         
!        IC:                 NUMBER OF RECORDED MODELS                  
!                                                                       
!xx      IMPLICIT  REAL*8  ( A-H,O-Z )                                     
!c      DIMENSION  X(500), F(501), PM(500), Z(500), C(20)                 
!c      DIMENSION  IX(500) , IND(500) , JND(1000) , KND(1000) , POST(1000)
!c      DIMENSION  Y(500)                                                 
!xx      DIMENSION  X(NN), F(NN+1), PM(NN), Z(NN), C(K+1)                 
!xx      DIMENSION  IX(NN) , IND(NN) , JND(2**K) , KND(2**K) , POST(2**K)
!xx      DIMENSION  Y(NN)                     
integer nn, k, isw, jsw, ioutd
real(dp) z(nn), y(nn), rlim, rout
! local
integer i, ii, ii1, ic, il, isprss, jj, jj1, k1, k2,&
&ix(nn), ind(nn), jnd(2**k), knd(2**k), n, nml1, l
real(dp) x(nn), f(nn+1), pm(nn), c(k+1), post(2**k), eps,&
&di, sumf, dlk0, f0, tem, ff, expf, w

!c      COMMON /CSPRSS/ ISPRSS                                            
isprss = 1
!                                                                       
n=0
do 5 i=1,nn
if(rlim .le. 0.d0) go to 6
if(z(i) .ge. rlim) go to 7
6 n=n+1
x(n)=z(i)
ix(n)=i
7 y(i)=z(i)
5 continue
l=k
!c      IF(ISPRSS .EQ. 0) WRITE( 6,600 )   N, K, L, ISW, JSW              
!c      IF(ISPRSS .EQ. 0) WRITE( 6,601 )   (X(I),I=1,N)                   
!                                                                       
!c      CALL  SRTMIN( X,N,IX )                                            
call  bsrtmin( x,n,ix )
eps = 1.0d-03
nml1 = n-l+1
f(1) = 0.d0
do 10  i=1,n
ind(i) = 0
pm(i) = 0.d0
di = i
!xx   10 F(I+1) = F(I)+DLOG(DI)                                            
f(i+1) = f(i)+dlog(di)
10 continue
!c      C(1) = DFLOAT(N*2)/DFLOAT(N-3)                                    
c(1) = dble(n*2)/dble(n-3)
do 20  i=1,k
!xx   20 C(I+1) = DFLOAT(N*(I+2))/DFLOAT(N-I-3)+F(N+1)-F(N-I+1)            
!c      C(I+1) = DFLOAT(N*(I+2))/DFLOAT(N-I-3)+F(N+1)-F(N-I+1)
c(i+1) = dble(n*(i+2))/dble(n-i-3)+f(n+1)-f(n-i+1)
20 continue
!                                                                       
il = 2**l
if(jsw.eq.1)  il=k+1
ic = 0
sumf = 0.d0
dlk0 = 0.d0
!                                                                       
do 101  ii1 = 1,il
ii = ii1-1
if(jsw.eq.1)  ii=2**ii-1
call  binary( ii,l,ind )
k1 = 0
do 30  i=1,l
!xx   30 K1 = K1+IND(I)                                                    
k1 = k1+ind(i)
30 continue
if( k1 .gt. k )   go to 101
!                                                                       
do 100  jj1=1,il
jj = jj1-1
if(jsw.eq.1)   jj=2**k-2**jj
call  binary( jj,l,ind(nml1) )
k2 = k1
do 40  i=nml1,n
!xx   40 K2 = K2+IND(I)                                                    
k2 = k2+ind(i)
40 continue
if( k2 .gt. k )   go to 100
!                                                                       
call  lkout1( x,n,ind,jsw,ff,w )
if( ic.eq.0 ) dlk0=ff-c(k2+1)
!                                                                       
f0 = ff-c(k2+1)-dlk0
if( f0 .lt. -20.d0 )  go to 100
if(f0.lt.20.d0)go to 45
dlk0=f0
tem=dexp(-f0)
sumf=sumf*tem
do 46 i=1,n
!xx   46 PM(I)=PM(I)*TEM                                                   
pm(i)=pm(i)*tem
46 continue
do 47 i=1,ic
!xx   47 POST(I)=POST(I)*TEM                                               
post(i)=post(i)*tem
47 continue
f0=0.d0
45 continue
expf = dexp( f0 )*w
sumf = sumf+expf
do 50  i=1,n
!xx   50 PM(I) = PM(I) + IND(I)*EXPF                                       
pm(i) = pm(i) + ind(i)*expf
50 continue
if( expf/sumf .lt. eps )   go to 100
!                                                                       
ic = ic+1
jnd(ic) = ii
knd(ic) = jj
post(ic) = expf
!                                                                       
100 continue
!                                                                       
101 continue
do 110  i=1,n
!xx  110 PM(I) = PM(I)/SUMF                                                
pm(i) = pm(i)/sumf
110 continue
!                                                                       
do 120  i=1,ic
!xx  120 POST(I) = POST(I)/SUMF                                            
post(i) = post(i)/sumf
120 continue
!                                                                       
call  prpost( post,x,ind,jnd,knd,ic,n,l)
!                                                                       
!c      IF(ISPRSS .EQ. 0) WRITE( 6,602 )   (PM(I),I=1,N)                  
!c      IF( ISW .GE. 1 )   CALL  MODIFY( N,L,IX,PM,JND,KND,Y,IC )         
if( isw .ge. 1 )  call modify( n,l,ix,pm,jnd,knd,y,ic,ioutd,rout )
!c      IF(ISW .GE. 1 .AND. ISPRSS .EQ. 0) WRITE(6,603)(Y(I),I=1,NN)      
!                                                                       
return
!xx  600 FORMAT( 1H ,'N    =',I6,5X,'(NUMBER OF DATA)',/,                  
!xx     1' K    =',I6,5X,'(MAXIMUM NUMBER OF OUTLIERS)',/,                 
!xx     2' L    =',I6,5X,'(RANGE OF SEARCH ON BOTH SIDES)',/,              
!xx     3' ISW  =',I6,/,' JSW  =',I6 )                                     
!xx  601 FORMAT( 1H ,'**  DATA  **',/,(1X,12D11.3) )                       
!xx  602 FORMAT( 1H ,'**  MARGINAL POSTERIOR PROBABILITIES  **',/,(1X,12D11
!xx     *.3) )                                                             
!xx  603 FORMAT(1H ,'** MODIFIED DATA **',/,(1X,12D11.3) )                 
end
!c      SUBROUTINE  MODIFY( N,L,IX,POST,JND,KND,Y,IC )                    
subroutine  modify( n,l,ix,post,jnd,knd,y,ic,ioutd,const )
  use timsac_kinds, only: dp
  implicit none
!     REVISED  MAY 29, 1980                                             
!                                                                       
!     THIS SUBROUTINE MODIFIES THE ORIGINAL DATA BY USING THE OUTPUTS   
!     FROM SUBROUTINE OUTLIR.                                           
!                                                                       
!     INPUTS:                                                           
!        (X(I),I=1,N):   ORIGINAL DATA                                  
!        N:              NUMBER OF ORIGINAL DATA                        
!        L:              NUMBER OF POSSIBLE OUTLIERS IN BOTH SIDES      
!        (IX(I),I=1,N):                                                 
!        (POST(I),I=1,N): MARGINAL POSTERIOR PROBABILITIES THAT X(IX(I))
!                        IS AN OUTLIER                                  
!        IOUTD:          =1  X(IX(I)) IS JUDGED AS AN OUTLIER           
!                            WHEN POST(I) IS GREATER THAN 0.01          
!                        =2  X(IX(I)) IS JUDGED AS AN OUTLIER           
!                            WHEN AT LEAST ONE MODEL, WHOSE             
!                            POSTERIOR PROBABILITY IS GREATER THAN      
!                            THAT OF NO-OUTLIER MODEL, SPECIFIES IT     
!                            AS AN OUTLIER                              
!                                                                       
!     OUTPUT:                                                           
!        (Y(I),I=1,N):   MODIFIED DATA                                  
!                                                                       
!xx      IMPLICIT  REAL*8  ( A-H,O-Z )                                     
!x      DIMENSION   IX(N),  POST(1), JND(1), KND(1), Y(N)                 
!c      DIMENSION  IND(500)                                               
!xx      DIMENSION   IX(N),  POST(N), JND(IC), KND(IC), Y(N)                 
!xx      DIMENSION  IND(N)              
integer n, l, ic ,ioutd, ix(n), jnd(ic), knd(ic)
real(dp) post(n), y(n), const
! local
integer i, ictem, ichk, j, k, nml1, ind(n)
!c      COMMON /CSPRSS/ ISPRSS                                            
!c      COMMON /ILOGT/ IDUMMY(3),IOUTD,CONST                              
ictem = ic
if(ioutd .eq. 1) ictem=1
!xx      DO 110 I=1,N                                                      
!xx  110 IND(I) = 0                                                        
ind(1:n) = 0
nml1 = n-l+1
do 200 k=1,ictem
call  binary( jnd(k),l,ind )
call  binary( knd(k),l,ind(nml1) )
!                                                                       
ichk = 0
do 120  i=1,n
j = ix(i)
if(ioutd .eq. 1 .and. post(i) .le. 0.01d0) go to 120
if(ioutd .eq. 2 .and. ind(i) .eq. 0) go to 120
ichk=1
y(j) = const
120 continue
if(ichk .eq. 0) go to 201
200 continue
201 continue
return
!xx  600 FORMAT( 1H ,'**  MODIFIED DATA (ISW=1)  **',/,(1X,10F13.5) )      
!xx  601 FORMAT( 1H ,'**  MODIFIED DATA (ISW=2)  **',/,(1X,10F13.5) )      
!xx  602 FORMAT( 1H ,'**  MODIFIED DATA (ISW=3)  **',/,(1X,10F13.5) )      
!xx  605 FORMAT( 1H ,'LOCATION PARAMETER;   XM =',F13.5 )                  
end
subroutine setd(w,ip,id,c,iar,ar)
  use timsac_kinds, only: dp
  implicit none
!xx      IMPLICIT REAL*8 (A-H,O-Z)                                         
!c      DIMENSION W(IP,10), WW(10), AR(10)                                
!xx      DIMENSION W(IP,ID+IAR+1), WW(ID+IAR+1), AR(1)
integer ip, id, iar, jlx
real(dp) w(ip,id+iar+1), c, ar(1), ww(id+iar+1)
! local
integer i, idar, idp1, j, jy, k
!
idar = id + iar
idp1 = idar + 1
w(1,idp1) = c
ww(idp1) = c
if(idar .eq. 0) go to 998
!xx      DO 10 J=1,IDAR                                                    
!xx      WW(J)=0.D0                                                        
!xx      DO 10 I=1,IP                                                      
!xx   10 W(I,J) = 0.D0
ww(1:idar)=0.d0
w(1:ip,1:idar) = 0.d0
if(id .eq. 0) go to 997
!xx      DO 50 J=1,ID                                                      
do 51 j=1,id
i = idp1 - j - 1
do 50 k=1,j
i = i + 1
!xx   50 WW(I) = WW(I) - WW(I+1)
ww(i) = ww(i) - ww(i+1)
50 continue
51 continue
997 continue
do 500 j=1,idar
w(1,j) = ww(j)
jlx = min0(iar,idp1 - j)
if(iar .eq. 0) go to 500
do 400 jy=1,jlx
!xx  400 W(1,J) = W(1,J) - AR(JY)*WW(J+JY)                                 
w(1,j) = w(1,j) - ar(jy)*ww(j+jy)
400 continue
500 continue
998 continue
return
end
subroutine  init(w,length,dop,istep)
  use timsac_kinds, only: dp
  implicit none
!xx      IMPLICIT REAL*8 (A-H,O-Z)                                         
!c      DIMENSION  W(100), DDOP(100), DOP(100)                            
!xx      DIMENSION  W(LENGTH), DDOP(LENGTH), DOP((LENGTH-1)*ISTEP+1)
integer length, istep
real(dp) w(length), dop((length-1)*istep+1)
! local
integer i, item, j, k
real(dp) ddop(length), sum
j=1
do 1 i=1,length
ddop(i)=dop(j)
!xx    1 J=J+ISTEP
j=j+istep
1 continue
do 20 j=1,length
sum = 0.d0
item = 0
do 10 k=j,length
item=item+1
!xx   10 SUM=SUM-W(K)*DDOP(ITEM)
sum=sum-w(k)*ddop(item)
10 continue
!xx   20 W(J)=SUM
w(j)=sum
20 continue
return
end
!                                                                       
!                                                                       
!                                                                       
subroutine  exhsld(h1,n1,h2,n2,h3,n3,h4,m1,ipos)
  use timsac_kinds, only: dp
  implicit none
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!x      DIMENSION  H1(N1,1),H2(N2,1),H3(1),H4(1)                          
!xx      DIMENSION  H1(N1,IPOS),H2(N2,N2+IPOS),H3(N3),H4(N2)                          
integer n1, n2, n3, m1, ipos
real(dp) h1(n1,ipos), h2(n2,n2+ipos), h3(n3), h4(n2)
! local
integer j, jp1, k, m, m0, mm
real(dp) eps, sqd, c, d, f
data eps/1.d-30/
if(ipos .le. m1) go to 30
m1 = ipos
!xx      DO 10 J=1,N1                                                      
!xx   10 H1(J,M1) = 0.D0
!xx      DO 20 J=1,N2                                                      
!xx   20 H2(J,M1+N2) = 0.D0
h1(1:n1,m1) = 0.d0
h2(1:n2,m1+n2) = 0.d0
30 continue
if(n3 .lt. 0) return
m0 = ipos - n3
do 100 j=1,n3
mm = m0 + j
if(dabs(h3(j)).lt.eps)go to 100
d = h1(1,mm)**2 + h3(j)**2
sqd = dsqrt(d)
if(h1(1,mm).gt.0.d0) sqd=-sqd
c = d - sqd*h1(1,mm)
d = h1(1,mm) - sqd
h1(1,mm) = sqd
jp1=j+1
if(jp1.gt.n3)go to 60
m = 1
do 50 k=jp1,n3
m = m + 1
if(m .gt. n1) go to 60
f = d*h1(m,mm) + h3(j)*h3(k)
f = f/c
h1(m,mm) = h1(m,mm) - d*f
!xx   50 H3(K) = H3(K) - H3(J)*F                                           
h3(k) = h3(k) - h3(j)*f
50 continue
60 continue
do 70 k=1,n2
f = d*h2(k,mm) + h3(j)*h4(k)
f = f/c
h2(k,mm) = h2(k,mm) - d*f
!xx   70 H4(K) = H4(K) - H3(J)*F                                           
h4(k) = h4(k) - h3(j)*f
70 continue
100 continue
do 200 j=1,n2
mm = m1 + j
if(dabs(h4(j)).lt.eps)go to 200
d = h2(j,mm)**2 + h4(j)**2
sqd = dsqrt(d)
if(h2(j,mm).gt.0.d0) sqd=-sqd
c = d - sqd*h2(j,mm)
d = h2(j,mm) - sqd
h2(j,mm) = sqd
if(j .ge. n2) go to 200
jp1 = j + 1
do 150 k=jp1,n2
f = d*h2(k,mm) + h4(j)*h4(k)
f = f/c
h2(k,mm) = h2(k,mm) - d*f
!xx  150 H4(K) = H4(K) - H4(J)*F                                           
h4(k) = h4(k) - h4(j)*f
150 continue
200 continue
return
end
!                                                                       
!                                                                       
!                                                                       
!c      SUBROUTINE SETX(H1,N1,H2,N2,M1,ICOUNT,FTRN,N,H,NH,WT,Y,NDATA,RLIM,
!c     *                 WEEK,IY,ALNDTD,F,DD,IART,ARFT)                   
subroutine setx(h1,n1,h2,n2,m1,icount,ftrn,n,h,nh,wt,y,&
&ndata,rlim,week,iy,alndtd,f,dd,iart,arft,&
!xx     *                   ALPHA,BETA,GAMMA,WTRD,DELTA,IP,ID,IS,YEAR,NPF)
&alpha,wtrd,delta,ip,id,year,npf)
  use timsac_kinds, only: dp
  implicit none
!xx      IMPLICIT REAL*8 (A-H,O-Z)                                         
!xx      INTEGER*4 YEAR                                                    
!c      DIMENSION H1(N1,1),H2(N2,1),Y(1),WEEK(IY,1),H3(200),H4(50)
!c      DIMENSION FTRN(1),H(NH,1),TI(10),T0(10),F(1),ARFT(1)              
!x      DIMENSION H1(N1,1),H2(N2,1),Y(NDATA),WEEK(IY,1),H3(N1),H4(N2)
!xx      DIMENSION H1(N1,1),H2(N2,N+8),Y(NDATA),WEEK(IY,N),H3(N1),H4(N2)
!xx      DIMENSION FTRN(ID+3),H(NH,N),TI(ID+IART),T0(2*(ID+IART+1)),
!xx     *           F(NPF+1),ARFT(3)
integer n1, n2, m1, icount, n, nh, ndata, iy, iart, ip, id,&
&year, npf
real(dp) h1(n1,1), h2(n2,n+8), ftrn(id+3), h(nh,n),&
&ti(id+iart), wt, y(ndata), rlim, week(iy,n),&
&alndtd, f(npf+1), dd, arft(3), alpha, wtrd, delta
! local
integer i, item, id0, idar, ipos, j, j0, jtem, k, kp1, n2m1,&
&n3, nmk
real(dp) h3(n1), h4(n2), t0(2*(id+iart+1)), tem, temo
!c      COMMON/IDATA/IP,ID,IS,YEAR                                        
!c      COMMON/ RDATA/ALPHA,BETA,GAMMA,DUMMY(4),WTRD,DELTA                
!xx      DO 600 I=1,N2                                                     
!xx  600 H2(I,N2)=0.D0
h2(1:n2,n2)=0.d0
if(year .eq. 0) go to 200
temo=wtrd
tem=-temo/7.d0
do 300 j=1,7
do 400 i=1,7
!xx  400 H2(I,J)=TEM                                                       
h2(i,j)=tem
400 continue
h2(8,j)=tem*delta
!xx  300 H2(J,J)=H2(J,J)+TEMO                                              
h2(j,j)=h2(j,j)+temo
300 continue
!c      CALL HUSHLD(H2,8,8,8,0)                                           
call bhushld(h2,8,8,8,0)
!xx      DO 500 I=1,N2                                                     
do 501 i=1,n2
do 500 j=1,i
tem=h2(i,j)
h2(i,j)=h2(j,i)
!xx  500 H2(J,I)=TEM                                                       
h2(j,i)=tem
500 continue
501 continue
n2m1=n2-1
do 550 j=1,n2m1
tem=dabs(h2(j,j))
!xx  550 ALNDTD=ALNDTD+DLOG(TEM)                                           
alndtd=alndtd+dlog(tem)
550 continue
200 continue
item=2
if(ip.eq.1)item=1
call setd(t0,item,id,wt,iart,arft)
idar = id + iart
do 10 i=1,idar
!xx   10 TI(I) = FTRN(I)*ALPHA
ti(i) = ftrn(i)*alpha
10 continue
call init(ti,idar,t0,item)
icount=0
m1=0
id0=idar*item+1
ipos=0
k=0
do 100 i=1,n
ipos = ipos + 1
n3 = ipos
if(n3 .gt. id0) n3 = id0
jtem = id0 - n3
do 110 j=1,n3
jtem = jtem + 1
!xx  110 H3(J) = T0(JTEM)                                                  
h3(j) = t0(jtem)
110 continue
!xx      DO 120 J=1,N2                                                     
!xx  120 H4(J) = 0.D0                                                      
h4(1:n2) = 0.d0
if(i .gt. idar) go to 125
h4(n2)=ti(i)
do 121 j=1,n3
!xx  121 H3(J)=H3(J)*ALPHA
h3(j)=h3(j)*alpha
121 continue
125 continue
call exhsld(h1,n1,h2,n2,h3,n3,h4,m1,ipos)
n3 = -1
if(ip .gt. 1) ipos=ipos+1
if(i .gt. ndata) go to 90
if(y(i) .gt. rlim .and. rlim .gt. 0.d0) go to 90
icount=icount+1
n3 = item
h3(1) = 1.d0
if(ip.ne.1)h3(2) = 1.d0
do 30 j=1,n2
if(j .lt. n2) h4(j) = week(j,i)
if(j .eq. n2) h4(j) = y(i)
30 continue
90 call exhsld(h1,n1,h2,n2,h3,n3,h4,m1,ipos)
if(ip .eq. 1) go to 100
if(ipos .le. n1) go to 100
k=k+1
j0=0
do 91 j=1,n1
if(mod(j,2) .eq. 0) go to 911
j0=j0+1
h3(j)=h(j0,k)*dd
go to 91
911 h3(j)=0.d0
91 continue
!xx      DO 92 J=1,N2                                                      
!xx   92 H4(J)=0.D0                                                        
h4(1:n2)=0.d0
h4(n2)=f(k)*dd
call exhsld(h1,n1,h2,n2,h3,n1,h4,m1,ipos)
100 continue
if(ip .eq. 1) return
nmk=n1
kp1=k+1
do 95 k=kp1,n
nmk=nmk-2
j0=0
do 96 j=1,nmk
if(mod(j,2).eq.0)go to 93
j0=j0+1
h3(j)=h(j0,k)*dd
go to 96
93 h3(j)=0.d0
96 continue
!xx      DO 97 L=1,N2                                                      
!xx   97 H4(L)=0.D0
h4(1:n2)=0.d0
h4(n2)=f(k)*dd
!xx   95 CALL EXHSLD(H1,N1,H2,N2,H3,NMK,H4,M1,M1)                          
call exhsld(h1,n1,h2,n2,h3,nmk,h4,m1,m1)
95 continue
return
end
!                                                                       
!                                                                       
!                                                                       
!c      SUBROUTINE SETDC(H1,N1,H2,M1,FSEAS,N,WS,WZ,IARS,ARFS,IARN,ARFN)   
!c     *                 ARFN,ALPHA,BETA,GAMMA,WTRD,DELTA,IP,ID,IS,YEAR)
subroutine setdc(h1,n1,h2,m1,fseas,n,ws,wz,iars,arfs,iarn,&
!x     *                  ARFN,ALPHA,BETA,GAMMA,WTRD,DELTA,IP,IS,YEAR)
!xx     *                  ARFN,ALPHA,BETA,GAMMA,WTRD,DELTA,IP,IS,YEAR,NPF)
&arfn,beta,gamma,ip,is,npf)
  use timsac_kinds, only: dp
  implicit none
!xx      IMPLICIT REAL*8 (A-H,O-Z)                                         
!xx      INTEGER*4 YEAR                                                    
!c      DIMENSION  H1(N1,1),H2(1),FSEAS(1),H3(100),S0(100),ARFS(1),ARFN(1)
!c      DIMENSION  Z0(50),SI(100),ZI(50),H4(50)                           
!x      DIMENSION  H1(N1,1),H2(1),FSEAS(1),H3((IS+IARS)*IP+IARN+1),
!x     *            S0((IS+IARS+1)*IP+IARN+1),ARFS(1),ARFN(1)
!x      DIMENSION  Z0(IP),SI((IS+IARS)*IP+IARN),ZI(IP),H4(1)
!xx      DIMENSION  H1(N1,1),H2(NPF+1),H3((IS+IARS)*IP+IARN+1),h4(1)
!xx      DIMENSION  FSEAS((IS+IARS)*IP + IARN),S0((IS+IARS+1)*IP+IARN+1)
!xx      DIMENSION  Z0(IP),SI((IS+IARS)*IP+IARN),ZI(IP),ARFS(3),ARFN(3)
integer n1, m1, n, iars, iarn, ip, is, npf
real(dp) h1(n1,1), h2(npf+1), fseas((is+iars)*ip + iarn),&
&ws, wz, arfs(3), arfn(3), beta, gamma
! local
integer i, ipis, ipm1, is0, ipos, item, iz0, j, ltem, length, n3
real(dp) h3((is+iars)*ip+iarn+1), h4(1), z0(ip),  zi(ip),&
&sum, s0((is+iars+1)*ip+iarn+1),&
&si((is+iars)*ip+iarn), dt
!c      COMMON /RDATA/ ALPHA,BETA,GAMMA,DUMMY(4),WTRD,DELTA               
!c      COMMON /IDATA/ IP,ID,IS,YEAR                                      
ipis = (is+iars)*ip + iarn
ipm1 =ip - 1
item = ip*(is-1) + 1
do 30 i=1,ipm1
item = item + 1
!xx   30 ZI(I) = FSEAS(ITEM)*WZ*GAMMA                                      
zi(i) = fseas(item)*wz*gamma
30 continue
sum = 0.d0
item = ip
do 40 i=1,ipm1
item = item - 1
sum = sum - zi(item)
!xx   40 ZI(ITEM) = SUM                                                    
zi(item) = sum
40 continue
call setd(s0,ip,is,ws,iars,arfs)
if(iarn .eq. 0) go to 49
length = ipis + 1
ltem = length - iarn
do 41 i=1,ltem
s0(length)=s0(length - iarn)
!xx   41 LENGTH = LENGTH - 1                                               
length = length - 1
41 continue
!xx      DO 42 I=1,IARN                                                    
!xx   42 S0(I)=0.D0
s0(1:iarn)=0.d0
length = ipis + 1
do 400 i=1,length
dt=s0(i)
do 410 j=1,iarn
!xx  410 IF(I+J .LE. LENGTH) DT=DT-ARFN(J)*S0(I+J)                         
if(i+j .le. length) dt=dt-arfn(j)*s0(i+j)
410 continue
!xx  400 S0(I)=DT                                                          
s0(i)=dt
400 continue
49 continue
if(ipis .eq. 0) go to 998
do 20 i=1,ipis
!xx   20 SI(I) = FSEAS(I)*BETA                                             
si(i) = fseas(i)*beta
20 continue
call init(si,ipis,s0,1)
998 do 50 i=1,ip
z0(i)=wz
50 continue
!xx   55 CONTINUE                                                          
is0 = ipis + 1
iz0 = ip
m1 = 0
ipos=0
h2(1) = 0.d0
do 100 i=1,n
ipos = i
n3 = ipos
if(n3 .gt. is0) n3 = is0
item = is0 - n3
do 130 j=1,n3
item = item + 1
!xx  130 H3(J) = S0(ITEM)                                                  
h3(j) = s0(item)
130 continue
h4(1) = 0.d0
if(i .gt. ipis) go to 145
h4(1)=si(i)
do 141 j=1,n3
!xx  141 H3(J)=H3(J)*BETA
h3(j)=h3(j)*beta
141 continue
145 continue
call exhsld(h1,n1,h2,1,h3,n3,h4,m1,ipos)
n3 = ipos
if(n3 .gt. iz0) n3 =iz0
item = iz0 - n3
do 150 j=1,n3
item = item + 1
!xx  150 H3(J) = Z0(ITEM)                                                  
h3(j) = z0(item)
150 continue
h4(1) = 0.d0
if(i .ge. ip) go to 165
h4(1)=zi(i)
do 161 j=1,n3
!xx  161 H3(J)=H3(J)*GAMMA                                                 
h3(j)=h3(j)*gamma
161 continue
165 continue
call exhsld(h1,n1,h2,1,h3,n3,h4,m1,ipos)
100 continue
return
end
subroutine partar(r,a,m)
  use timsac_kinds, only: dp
  implicit none
! **** PARCOR R TO AR ********                                          
!xx      IMPLICIT REAL*8(A-H,O-Z)                                          
!xx      DIMENSION R(M),A(M,M)                                             
integer m
real(dp) r(m), a(m,m)
! local
integer i, im1, j
!xx      DO 10 I=1,M                                                       
do 20 i=1,m
do 10 j=1,i
!xx   10 A(I,J)=0.D0
a(i,j)=0.d0
10 continue
20 continue
a(1,1)=r(1)
if(m.le.1) return
do 40 i=2,m
a(i,i)=r(i)
im1=i-1
do 30 j=1,im1
!xx   30 A(I,J)=A(I-1,J)-R(I)*A(I-1,I-J)                                   
a(i,j)=a(i-1,j)-r(i)*a(i-1,i-j)
30 continue
40 continue
return
end
!c      SUBROUTINE SOLVE(H1,N1,H2,N2,A,M1,SQE,NANS,ERR)                   
subroutine bsolve(h1,n1,h2,n2,a,m1,sqe,nans,err)
  use timsac_kinds, only: dp
  implicit none
!xx      IMPLICIT REAL*8 (A-H,O-Z)                                         
!x      DIMENSION  H1(N1,1), H2(N2,1), A(1), ERR(1)                       
!xx      DIMENSION  H1(N1,N2+M1), H2(N2,N2+M1), A(NANS), ERR(NANS)
integer n1, n2, m1 ,nans
real(dp) h1(n1,n2+m1), h2(n2,n2+m1), a(nans), sqe,&
&err(nans)
! local
integer i, j, jj, k, kk, ka, kam1, kka, km1, kcopy, katem,&
&l, ler, ltem, n2mj
real(dp) aka
!                                                                       
!xx      DO 30 I=1,NANS                                                    
!xx   30 ERR(I)=0.D0
err(1:nans)=0.d0
do 100 ler = 1, nans
k = m1 + n2
ka = nans
jj = nans - 1
if( ler.eq.nans )  go to 44
sqe = 0.d0
kam1 = ka - 1
!xx            DO 40 KKA = 1, KAM1                                         
!xx   40       A(KKA) = 0.D0
a(1:kam1) = 0.d0
a(ler) = 1.d0
go to 48
44 continue
sqe = h2(n2,k)**2
kka = ka
kk  = k
do 46 j = 1, jj
kk = kk - 1
kka = kka - 1
a(kka) = h2(n2,kk)
46 continue
48 continue
!                                                                       
kcopy = k
do 50 j = 1, jj
ka = ka - 1
if( a(ka).eq.0.d0 )  go to 50
k = kcopy - j
if( j.ge.n2 )  go to 20
!                                                                       
a(ka) = a(ka)/h2(n2-j,k)
if( ler.lt.nans )  err(ka) = err(ka) + a(ka)**2
km1 = ka - 1
if( km1.le.0 )  go to 50
ltem = k
katem = ka
aka = a(ka)
n2mj = n2-j
do 10 l = 1, km1
katem = katem - 1
ltem = ltem - 1
a(katem) = a(katem) - aka*h2(n2mj,ltem)
10 continue
go to 50
!                                                                       
20 a(ka) = a(ka)/h1(1,k)
if( ler.lt.nans )  err(ka) = err(ka) + a(ka)**2
ltem = k
l = ka
if( n1.lt.2 )  go to 50
do 25 i = 2, n1
l = l - 1
ltem = ltem - 1
if( l.le.0 )  go to 50
a(l) = a(l) - a(ka)*h1(i,ltem)
25 continue
50 continue
100 continue
!                                                                       
return
end
