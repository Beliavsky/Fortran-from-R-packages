# SPDX-License-Identifier: GPL-2.0-or-later
# Genuine upstream TTR 0.24.4 oracle; run from the package root with Rscript test/reference.R.
suppressPackageStartupMessages(library(TTR))
stopifnot(as.character(packageVersion('TTR')) == '0.24.4')
N <- 96L
i <- seq_len(N)
x <- 100 + (i %% 17) * .75 + i * .125 + ((i %% 5) - 2) * .2
y <- 80 + (i %% 13) * .6 + i * .2
v <- 1000 + i * 17 + (i %% 7) * 100
hlc <- cbind(x + 2 + (i %% 3) * .1, x - 1 - (i %% 4) * .1, x)
hl <- hlc[,1:2]
one <- matrix(x, ncol=1)
xna <- c(NA_real_, NA_real_, x[-c(1,2)])
ties <- rep(c(1, 2, 2, 4, 3, 2, 1, 1), length.out=N)
cases <- list()
add <- function(name, value, expression) {
  if (!is.numeric(value)) stop(name)
  cases[[length(cases)+1]] <<- list(name=name, value=as.numeric(value), expr=expression)
}
add('SMA', SMA(x,5), 'sma(x, 5)')
add('EMA', EMA(x,5), 'ema(x, 5)')
add('EMA Wilder', EMA(xna,5,wilder=TRUE), 'ema(xna, 5, wilder = .true.)')
add('EMA ratio', suppressWarnings(EMA(x,5,ratio=.3)), 'ema(x, 5, ratio = 0.3_dp)')
add('DEMA', DEMA(x,5), 'dema(x, 5)')
add('DEMA factor', DEMA(xna,5,v=.4,wilder=TRUE), 'dema(xna, 5, v = 0.4_dp, wilder = .true.)')
add('WMA', WMA(x,5), 'wma(x, 5)')
add('WMA fixed', WMA(xna,5,c(1,3,2,5,4)), 'wma(xna, 5, [1.0_dp, 3.0_dp, 2.0_dp, 5.0_dp, 4.0_dp])')
add('VWMA', VWMA(x,v,5), 'vwma(x, v, 5)')
add('VWAP', VWAP(x,v,5), 'vwap(x, v, 5)')
add('EVWMA', EVWMA(x,v,5), 'evwma(x, v, 5)')
add('EVWMA constant', EVWMA(x,rep(1000,N),5), 'evwma(x, [1000.0_dp], 5)')
add('ZLEMA odd', ZLEMA(x,5), 'zlema(x, 5)')
add('ZLEMA even', ZLEMA(xna,6), 'zlema(xna, 6)')
add('HMA', HMA(x,9), 'hma(x, 9)')
add('ALMA', ALMA(x,9), 'alma(x, 9)')
add('ALMA parameters', ALMA(xna,6,offset=.4,sigma=3), 'alma(xna, 6, offset = 0.4_dp, sigma = 3.0_dp)')
for (name in c('Sum','Mean','Min','Max','Median','SD','Var')) {
  fn <- get(paste0('run',name)); fname <- paste0('run_',tolower(name))
  add(paste0('run',name), fn(xna,n=5), sprintf('%s(xna, 5)',fname))
  add(paste0('run',name,' cumulative'), fn(xna,n=5,cumulative=TRUE), sprintf('%s(xna, 5, cumulative = .true.)',fname))
}
add('runMedian upper', runMedian(ties,4,non.unique='max'), 'run_median(ties, 4, non_unique = 1)')
add('runMedian lower', runMedian(ties,4,non.unique='min'), 'run_median(ties, 4, non_unique = -1)')
add('runCov', runCov(xna,y,5), 'run_cov(xna, y, 5)')
add('runCov expanding population', runCov(xna,y,5,sample=FALSE,cumulative=TRUE),
    'run_cov(xna, y, 5, sample = .false., cumulative = .true.)')
add('runCor', runCor(x,y,5), 'run_cor(x, y, 5)')
add('runVar cross', runVar(x,y,n=5), 'run_var(x, 5, y = y)')
add('runSD population', runSD(x,5,sample=FALSE), 'run_sd(x, 5, sample = .false.)')
add('runMAD', runMAD(xna,5), 'run_mad(xna, 5)')
add('runMAD cumulative', runMAD(ties,4,cumulative=TRUE,non.unique='max'),
    'run_mad(ties, 4, cumulative = .true., non_unique = 1)')
add('runMAD mean', runMAD(x,5,stat='mean',center=SMA(x,5)),
    'run_mad(x, 5, center = sma(x, 5), median = .false.)')
add('runPercentRank ties', runPercentRank(ties,5), 'run_percent_rank(ties, 5)')
add('runPercentRank cumulative', runPercentRank(xna,5,cumulative=TRUE,exact.multiplier=1),
    'run_percent_rank(xna, 5, cumulative = .true., exact_multiplier = 1.0_dp)')
add('wilderSum', wilderSum(xna,5), 'wilder_sum(xna, 5)')
add('ROC continuous', ROC(x,3), 'roc(x, 3)')
add('ROC discrete', ROC(x,3,type='discrete'), 'roc(x, 3, .true.)')
add('momentum', momentum(x,3), 'momentum(x, 3)')
add('lags', lags(hlc,2), 'reshape(lags(hlc, 2), [94 * 9])')
add('growth default', growth(x), 'growth(x)')
add('growth discrete unpadded', growth(x,signals=rep(1,N-1),type='discrete',na.pad=FALSE), 'growth(x, discrete = .true., na_pad = .false.)')
for (name in c('TR','CLV','williamsAD')) {
 fn <- get(name); ff <- c(TR='tr',CLV='clv',williamsAD='williams_ad')[[name]]
 expr <- sprintf('%s(hlc)',ff)
 if(name=='TR') expr <- paste0('reshape(',expr,', [96 * 3])')
 add(name,fn(hlc),expr)
}
for (name in c('ATR','ADX','BBands','CCI','VHF','WPR','SNR')) {
 fn <- get(name); expr <- sprintf('%s(hlc, 5)',tolower(name)); nr <- c(ATR=4,ADX=4,BBands=4)
 if (name %in% names(nr)) expr <- paste0('reshape(',expr,', [96 * ',nr[[name]],'])')
 add(name,fn(hlc,5),expr)
}
for (name in c('CMO','RSI')) add(name,get(name)(x,5),sprintf('%s(x, 5)',tolower(name)))
add('CMF', CMF(hlc,v,5), 'cmf(hlc, v, 5)')
add('DPO', DPO(x,5), 'dpo(x, 5)')
add('DPO percent', DPO(x,6,shift=2,percent=TRUE,maType='EMA'),
    "dpo(x, 6, shift_by = 2, percent = .true., ma_type = 'EMA')")
add('DVI', DVI(x,8,smooth=3,magnitude=c(3,8,3),stretch=c(3,8,2)),
    'reshape(dvi(x, 8, smooth = 3, magnitude = [3, 8, 3], stretch = [3, 8, 2]), [96 * 3])')
add('DonchianChannel', DonchianChannel(hl,5), 'reshape(donchian_channel(hl, 5), [96 * 3])')
add('DonchianChannel lag', DonchianChannel(one,5,include.lag=TRUE),
    'reshape(donchian_channel(one, 5, include_lag = .true.), [96 * 3])')
add('EMV', EMV(hl,v,5), 'reshape(emv(hl, v, 5), [96 * 2])')
add('GMMA', GMMA(x,short=c(3,5),long=c(8,12)), 'reshape(gmma(x, [3, 5, 8, 12]), [96 * 4])')
add('KST', KST(x,n=c(3,4,5),nROC=c(2,4,6),nSig=3,wts=c(1,2,3)),
    'reshape(kst(x, [3, 4, 5], [2, 4, 6], 3), [96 * 2])')
add('MACD', MACD(x,3,8,4), 'reshape(macd(x, 3, 8, 4), [96 * 2])')
add('MACD difference SMA', MACD(x,3,8,4,percent=FALSE,maType='SMA'),
    "reshape(macd(x, 3, 8, 4, percent = .false., ma_type = 'SMA'), [96 * 2])")
add('MFI', MFI(hlc,v,5), 'mfi(hlc, v, 5)')
add('OBV', OBV(x,v), 'obv(x, v)')
add('OBV small changes', OBV(100+i*1e-8,v), 'obv(100.0_dp + [(real(i, dp) * 1.0e-8_dp, i = 1, 96)], v)')
add('OBV missing prefix', OBV(xna,c(NA,NA,v[-c(1,2)])), 'obv(xna, [nan_value(), nan_value(), v(3:)])')
add('PBands', PBands(x,5), 'reshape(pbands(x, 5), [96 * 3])')
add('PBands centered', PBands(x,5,centered=TRUE,lavg=TRUE),
    'reshape(pbands(x, 5, centered = .true., long_average = .true.), [96 * 3])')
add('TDI', TDI(x,5), 'reshape(tdi(x, 5), [96 * 2])')
add('TRIX', TRIX(x,5,3), 'reshape(trix(x, 5, 3), [96 * 2])')
add('chaikinAD', chaikinAD(hlc,v), 'chaikin_ad(hlc, v)')
add('chaikinVolatility', chaikinVolatility(hl,5), 'chaikin_volatility(hl, 5)')
add('keltnerChannels', keltnerChannels(hlc,5), 'reshape(keltner_channels(hlc, 5), [96 * 3])')
add('stoch', stoch(hlc,5,3,2), 'reshape(stoch(hlc, 5, 3, 2), [96 * 3])')
add('stoch unbounded smoothed', stoch(hlc,5,3,2,bounded=FALSE,smooth=2),
    'reshape(stoch(hlc, 5, 3, 2, bounded = .false., smooth = 2), [96 * 3])')
add('SMI', SMI(hlc,5,2,4,3), 'reshape(smi(hlc, 5, 2, 4, 3), [96 * 2])')
add('SMI unbounded', SMI(hlc,5,2,4,3,bounded=FALSE), 'reshape(smi(hlc, 5, 2, 4, 3, bounded = .false.), [96 * 2])')
add('ultimateOscillator', ultimateOscillator(hlc,c(3,5,8)),
    'ultimate_oscillator(hlc, [3, 5, 8])')
# Missing-value validation is computational, even though R returns the result invisibly.
nac <- naCheck(xna,5)
stopifnot(nac$NAs==2, nac$beg==7, identical(nac$nonNA,3:96))
con <- file('test/reference.dat','w')
for (cs in cases) {
 writeLines(cs$name,con)
 writeLines(as.character(length(cs$value)),con)
 writeLines(ifelse(is.na(cs$value),'NaN',sprintf('%.17g',cs$value)),con)
}
close(con)
header <- c('! SPDX-License-Identifier: GPL-2.0-or-later',
 '! Generated by test/reference.R against genuine TTR 0.24.4.',
 'program test_reference', '   use ttr_base', '   use ttr_indicators',
 '   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite',
 '   implicit none',
 '   real(dp) :: x(96), y(96), v(96), hlc(96, 3), hl(96, 2), one(96, 1), xna(96), ties(96)',
 '   integer :: i, unit, missing, begin', '   logical :: valid',
 '   do i = 1, 96',
 '      x(i) = 100.0_dp + real(mod(i, 17), dp) * 0.75_dp + real(i, dp) * 0.125_dp &',
 '         + real(mod(i, 5) - 2, dp) * 0.2_dp',
 '      y(i) = 80.0_dp + real(mod(i, 13), dp) * 0.6_dp + real(i, dp) * 0.2_dp',
 '      v(i) = 1000.0_dp + real(i * 17 + mod(i, 7) * 100, dp)',
 '      hlc(i, 1) = x(i) + 2.0_dp + real(mod(i, 3), dp) * 0.1_dp',
 '      hlc(i, 2) = x(i) - 1.0_dp - real(mod(i, 4), dp) * 0.1_dp',
 '      hlc(i, 3) = x(i)', '   end do',
 '   hl = hlc(:, :2)', '   one(:, 1) = x', '   xna = x', '   xna(:2) = nan_value()',
 '   do i = 1, 96, 8', '      ties(i:i + 7) = [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, 1.0_dp]', '   end do',
 '   open(newunit = unit, file = "test/reference.dat", status = "old", action = "read")')
body <- unlist(lapply(cases,function(cs)paste0('   call check(',cs$expr,')')))
footer <- c('   close(unit)', '   call na_check(xna, 5, missing, begin, valid)',
 '   if (missing /= 2 .or. begin /= 7 .or. .not. valid) error stop "naCheck fixture"',
 sprintf('   print *, "TTR genuine R oracle: %d numerical cases plus naCheck passed"',length(cases)),
 'contains', '   subroutine check(actual)',
 '      real(dp), intent(in) :: actual(:) !! Flattened translated output, compared with the next R fixture.',
 '      real(dp), allocatable :: expected(:)',
 '      integer :: m, j', '      character(120) :: label',
 '      read(unit, "(a)") label', '      read(unit, *) m',
 '      if (size(actual) /= m) error stop "reference shape"', '      allocate(expected(m))',
 '      do j = 1, m', '         read(unit, *) expected(j)', '      end do',
 '      do j = 1, m',
 '         if (ieee_is_nan(expected(j))) then',
 '            if (ieee_is_nan(actual(j))) cycle',
 '         else if (ieee_is_finite(expected(j)) .and. ieee_is_finite(actual(j))) then',
 '            if (abs(actual(j) - expected(j)) <= 2.0e-9_dp * max(1.0_dp, abs(expected(j)))) cycle',
 '         else if (.not. ieee_is_nan(actual(j))) then',
 '            if (actual(j) == expected(j)) cycle',
 '         end if',
 '         print *, trim(label), " row ", j, " actual ", actual(j), " expected ", expected(j)',
 '         error stop "R reference mismatch"', '      end do', '   end subroutine', 'end program')
writeLines(c(header,body,footer),'test/test_reference.f90')
cat(sprintf('Generated %d cases, %d expected numerical values.\n',length(cases),sum(vapply(cases,function(z)length(z$value),0L))))
