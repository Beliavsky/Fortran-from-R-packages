! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
program test_reference
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use vares
  implicit none
  integer :: failures
  failures=0
  call check_close('varexponential', varexponential(0.37_dp, lambda=1.3_dp), 0.35541189199735279_dp, 8.0e-08_dp, &
    & failures)
  call check_close('dexponential', dexponential(0.35541189199735279_dp, lambda=1.3_dp), 0.81900000000000006_dp, &
    & 2.0e-08_dp, failures)
  call check_close('pexponential', pexponential(0.35541189199735279_dp, lambda=1.3_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varkumexp', varkumexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.1398343293855876_dp, 8.0e-08_dp, failures)
  call check_close('dkumexp', dkumexp(0.1398343293855876_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 1.7449120092909294_dp, 2.0e-08_dp, failures)
  call check_close('pkumexp', pkumexp(0.1398343293855876_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varexpexp', varexpexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.26190704707873713_dp, 8.0e-08_dp, failures)
  call check_close('dexpexp', dexpexp(0.26190704707873713_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.94866895228000869_dp, 2.0e-08_dp, failures)
  call check_close('pexpexp', pexpexp(0.26190704707873713_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varbetaexp', varbetaexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.14921057184515724_dp, 8.0e-08_dp, failures)
  call check_close('dbetaexp', dbetaexp(0.14921057184515724_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 1.6621803573351386_dp, 2.0e-08_dp, failures)
  call check_close('pbetaexp', pbetaexp(0.14921057184515724_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varlogisexp', varlogisexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.31911059261090757_dp, 8.0e-08_dp, failures)
  call check_close('dlogisexp', dlogisexp(0.31911059261090757_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.71394294256559743_dp, 2.0e-08_dp, failures)
  call check_close('plogisexp', plogisexp(0.31911059261090757_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varexpext', varexpext(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.46743981345297458_dp, 8.0e-08_dp, failures)
  call check_close('dexpext', dexpext(0.46743981345297458_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.5958465269203429_dp, 2.0e-08_dp, failures)
  call check_close('pexpext', pexpext(0.46743981345297458_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.37000000000000011_dp, 2.0e-08_dp, failures)
  call check_close('varmoexp', varmoexp(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp), 0.44669768408773625_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dmoexp', dmoexp(0.44669768408773625_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.92219400000000018_dp, 2.0e-08_dp, failures)
  call check_close('pmoexp', pmoexp(0.44669768408773625_dp, lambda=1.3_dp, a=0.80000000000000004_dp), &
    & 0.36999999999999988_dp, 2.0e-08_dp, failures)
  call check_close('varperks', varperks(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.49540162336496069_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dperks', dperks(0.49540162336496069_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.40950000000000009_dp, 2.0e-08_dp, failures)
  call check_close('pperks', pperks(0.49540162336496069_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varbeard', varbeard(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, rho=0.59999999999999998_dp), &
    & 0.47904657033972764_dp, 8.0e-08_dp, failures)
  call check_close('dbeard', dbeard(0.47904657033972764_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & rho=0.59999999999999998_dp), 0.54609625782261495_dp, 2.0e-08_dp, failures)
  call check_close('pbeard', pbeard(0.47904657033972764_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & rho=0.59999999999999998_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('vargompertz', vargompertz(0.37_dp, b=1.7_dp, eta=1.1000000000000001_dp), &
    & 0.20628210162455526_dp, 8.0e-08_dp, failures)
  call check_close('dgompertz', dgompertz(0.20628210162455526_dp, b=1.7_dp, eta=1.1000000000000001_dp), &
    & 1.6729399772279141_dp, 2.0e-08_dp, failures)
  call check_close('pgompertz', pgompertz(0.20628210162455526_dp, b=1.7_dp, eta=1.1000000000000001_dp), &
    & 0.37000000000000011_dp, 2.0e-08_dp, failures)
  call check_close('varbetagompertz', varbetagompertz(0.37_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp, &
    & eta=1.1000000000000001_dp), 0.17882048111610946_dp, 8.0e-08_dp, failures)
  call check_close('dbetagompertz', dbetagompertz(0.17882048111610946_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp, &
    & eta=1.1000000000000001_dp), 1.1252724584238916_dp, 2.0e-08_dp, failures)
  call check_close('pbetagompertz', pbetagompertz(0.17882048111610946_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp, &
    & eta=1.1000000000000001_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varlfr', varlfr(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.40406850659629934_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dlfr', dlfr(0.40406850659629934_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.93675737056463682_dp, 2.0e-08_dp, failures)
  call check_close('plfr', plfr(0.40406850659629934_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.36999999999999988_dp, 2.0e-08_dp, failures)
  call check_close('varpareto', varpareto(0.37_dp, k=1.0_dp, c=1.2_dp), 1.469657748691884_dp, 8.0e-08_dp, failures)
  call check_close('dpareto', dpareto(1.469657748691884_dp, k=1.0_dp, c=1.2_dp), 0.51440548023708377_dp, &
    & 2.0e-08_dp, failures)
  call check_close('ppareto', ppareto(1.469657748691884_dp, k=1.0_dp, c=1.2_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varkumpareto', varkumpareto(0.37_dp, k=1.0_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), &
    & 1.1635633966622509_dp, 8.0e-08_dp, failures)
  call check_close('dkumpareto', dkumpareto(1.1635633966622509_dp, k=1.0_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & c=1.2_dp), 1.3842718095091289_dp, 2.0e-08_dp, failures)
  call check_close('pkumpareto', pkumpareto(1.1635633966622509_dp, k=1.0_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & c=1.2_dp), 0.36999999999999977_dp, 2.0e-08_dp, failures)
  call check_close('varF', varf(0.37_dp, d1=5.0_dp, d2=8.0_dp), 0.71456877568271637_dp, 8.0e-08_dp, failures)
  call check_close('dF', df(0.71456877568271637_dp, d1=5.0_dp, d2=8.0_dp), 0.61083287648167095_dp, 2.0e-08_dp, &
    & failures)
  call check_close('pF', pf(0.71456877568271637_dp, d1=5.0_dp, d2=8.0_dp), 0.36999999999999977_dp, 2.0e-08_dp, &
    & failures)
  call check_close('vargenpareto', vargenpareto(0.37_dp, k=0.40000000000000002_dp, c=1.2_dp), &
    & 0.14186921129302205_dp, 8.0e-08_dp, failures)
  call check_close('dgenpareto', dgenpareto(0.14186921129302205_dp, k=0.40000000000000002_dp, c=1.2_dp), &
    & 2.7420280854993044_dp, 2.0e-08_dp, failures)
  call check_close('pgenpareto', pgenpareto(0.14186921129302205_dp, k=0.40000000000000002_dp, c=1.2_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varbetapareto', varbetapareto(0.37_dp, k=1.0_dp, a=0.80000000000000004_dp, c=1.2_dp, d=1.5_dp), &
    & 1.6298456145391376_dp, 8.0e-08_dp, failures)
  call check_close('dbetapareto', dbetapareto(1.6298456145391376_dp, k=1.0_dp, a=0.80000000000000004_dp, c=1.2_dp, &
    & d=1.5_dp), 0.41372241611605026_dp, 2.0e-08_dp, failures)
  call check_close('pbetapareto', pbetapareto(1.6298456145391376_dp, k=1.0_dp, a=0.80000000000000004_dp, c=1.2_dp, &
    & d=1.5_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varparetostable', varparetostable(0.37_dp, lambda=1.3_dp, nu=4.5_dp, sigma=1.3_dp), &
    & 2.877696234177435_dp, 8.0e-08_dp, failures)
  call check_close('dparetostable', dparetostable(2.877696234177435_dp, lambda=1.3_dp, nu=4.5_dp, sigma=1.3_dp), &
    & 0.57282341947575066_dp, 2.0e-08_dp, failures)
  call check_close('pparetostable', pparetostable(2.877696234177435_dp, lambda=1.3_dp, nu=4.5_dp, sigma=1.3_dp), &
    & 0.37000000000000011_dp, 2.0e-08_dp, failures)
  call check_close('vargamma', vargamma(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.18376793322094942_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dGamma', dgamma(0.18376793322094942_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 1.348312412904159_dp, 2.0e-08_dp, failures)
  call check_close('pGamma', pgamma(0.18376793322094942_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varkumgamma', varkumgamma(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp), &
    & 0.15580737263330821_dp, 8.0e-08_dp, failures)
  call check_close('dkumgamma', dkumgamma(0.15580737263330821_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), 1.8074170769378839_dp, 2.0e-08_dp, failures)
  call check_close('pkumgamma', pkumgamma(0.15580737263330821_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varrgamma', varrgamma(0.37_dp, a=0.80000000000000004_dp, theta=0.59999999999999998_dp, &
    & phi=1.2_dp), 0.37411872289910647_dp, 8.0e-08_dp, failures)
  call check_close('drgamma', drgamma(0.37411872289910647_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp, phi=1.2_dp), 0.41405794386721961_dp, 2.0e-08_dp, failures)
  call check_close('prgamma', prgamma(0.37411872289910647_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp, phi=1.2_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varloggamma', varloggamma(0.37_dp, a=0.80000000000000004_dp, r=1.3999999999999999_dp), &
    & 0.16174209597900227_dp, 8.0e-08_dp, failures)
  call check_close('dloggamma', dloggamma(0.16174209597900227_dp, a=0.80000000000000004_dp, &
    & r=1.3999999999999999_dp), 1.509052396185566_dp, 2.0e-08_dp, failures)
  call check_close('ploggamma', ploggamma(0.16174209597900227_dp, a=0.80000000000000004_dp, &
    & r=1.3999999999999999_dp), 0.36999999999999977_dp, 2.0e-08_dp, failures)
  call check_close('varinvgamma', varinvgamma(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 2.2401606780467924_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dinvgamma', dinvgamma(2.2401606780467924_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.32249447017914806_dp, 2.0e-08_dp, failures)
  call check_close('pinvgamma', pinvgamma(2.2401606780467924_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varbetadist', varbetadist(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.17632046899030579_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dbetadist', dbetadist(0.17632046899030579_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 1.552303082372136_dp, 2.0e-08_dp, failures)
  call check_close('pbetadist', pbetadist(0.17632046899030579_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varuniform', varuniform(0.37_dp, a=0.0_dp, b=2.0_dp), 0.73999999999999999_dp, 8.0e-08_dp, &
    & failures)
  call check_close('duniform', duniform(0.73999999999999999_dp, a=0.0_dp, b=2.0_dp), 0.5_dp, 2.0e-08_dp, failures)
  call check_close('puniform', puniform(0.73999999999999999_dp, a=0.0_dp, b=2.0_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('vargenunif', vargenunif(0.37_dp, a=0.0_dp, c=1.0_dp, h=0.80000000000000004_dp, k=1.2_dp), &
    & 0.36560417929902111_dp, 8.0e-08_dp, failures)
  call check_close('dgenunif', dgenunif(0.36560417929902111_dp, a=0.0_dp, c=1.0_dp, h=0.80000000000000004_dp, &
    & k=1.2_dp), 1.0775466862667_dp, 2.0e-08_dp, failures)
  call check_close('pgenunif', pgenunif(0.36560417929902111_dp, a=0.0_dp, c=1.0_dp, h=0.80000000000000004_dp, &
    & k=1.2_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varpower1', varpower1(0.37_dp, a=0.80000000000000004_dp), 0.2885706482644807_dp, 8.0e-08_dp, &
    & failures)
  call check_close('dpower1', dpower1(0.2885706482644807_dp, a=0.80000000000000004_dp), 1.0257453479076992_dp, &
    & 2.0e-08_dp, failures)
  call check_close('ppower1', ppower1(0.2885706482644807_dp, a=0.80000000000000004_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varpower2', varpower2(0.37_dp, b=1.7_dp), 0.23798235276953797_dp, 8.0e-08_dp, failures)
  call check_close('dpower2', dpower2(0.23798235276953797_dp, b=1.7_dp), 1.4054792613957536_dp, 2.0e-08_dp, failures)
  call check_close('ppower2', ppower2(0.23798235276953797_dp, b=1.7_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varlogbeta', varlogbeta(0.37_dp, a=0.5_dp, b=1.5_dp, c=1.2_dp, d=2.0_dp), &
    & 1.2545227883117034_dp, 8.0e-08_dp, failures)
  call check_close('dlogbeta', dlogbeta(1.2545227883117034_dp, a=0.5_dp, b=1.5_dp, c=1.2_dp, d=2.0_dp), &
    & 3.2184605628032745_dp, 2.0e-08_dp, failures)
  call check_close('plogbeta', plogbeta(1.2545227883117034_dp, a=0.5_dp, b=1.5_dp, c=1.2_dp, d=2.0_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varcompbeta', varcompbeta(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.62410234344238946_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dcompbeta', dcompbeta(0.62410234344238946_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.90135461231237457_dp, 2.0e-08_dp, failures)
  call check_close('pcompbeta', pcompbeta(0.62410234344238946_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.37000000000000027_dp, 2.0e-08_dp, failures)
  call check_close('varLNbeta', varlnbeta(0.37_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.141383949440156_dp, 8.0e-08_dp, failures)
  call check_close('dLNbeta', dlnbeta(0.141383949440156_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 1.8571133584800945_dp, 2.0e-08_dp, failures)
  call check_close('pLNbeta', plnbeta(0.141383949440156_dp, lambda=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.37000000000000038_dp, 2.0e-08_dp, failures)
  call check_close('varMRbeta', varmrbeta(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, r=1.3999999999999999_dp, &
    & q=1.3_dp), 0.54978708485032723_dp, 8.0e-08_dp, failures)
  call check_close('dMRbeta', dmrbeta(0.54978708485032723_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & r=1.3999999999999999_dp, q=1.3_dp), 0.69696786457406579_dp, 2.0e-08_dp, failures)
  call check_close('pMRbeta', pmrbeta(0.54978708485032723_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & r=1.3999999999999999_dp, q=1.3_dp), 0.37000000000000038_dp, 2.0e-08_dp, failures)
  call check_close('vargenbeta', vargenbeta(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp), &
    & 1.2528961406970918_dp, 8.0e-08_dp, failures)
  call check_close('dgenbeta', dgenbeta(1.2528961406970918_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), 5.174343607907117_dp, 2.0e-08_dp, failures)
  call check_close('pgenbeta', pgenbeta(1.2528961406970918_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), 0.37000000000000044_dp, 2.0e-08_dp, failures)
  call check_close('vararcsine', vararcsine(0.37_dp, a=0.0_dp, b=2.0_dp), 0.60285210936521949_dp, 8.0e-08_dp, &
    & failures)
  call check_close('darcsine', darcsine(0.60285210936521949_dp, a=0.0_dp, b=2.0_dp), 0.34683550186038198_dp, &
    & 2.0e-08_dp, failures)
  call check_close('parcsine', parcsine(0.60285210936521949_dp, a=0.0_dp, b=2.0_dp), 0.37000000000000005_dp, &
    & 2.0e-08_dp, failures)
  call check_close('vartriangular', vartriangular(0.37_dp, a=0.0_dp, b=2.0_dp, c=0.80000000000000004_dp), &
    & 0.76941536246685382_dp, 8.0e-08_dp, failures)
  call check_close('dtriangular', dtriangular(0.76941536246685382_dp, a=0.0_dp, b=2.0_dp, &
    & c=0.80000000000000004_dp), 0.96176920308356728_dp, 2.0e-08_dp, failures)
  call check_close('ptriangular', ptriangular(0.76941536246685382_dp, a=0.0_dp, b=2.0_dp, &
    & c=0.80000000000000004_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('vargenbeta2', vargenbeta2(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), &
    & 0.23546097803593721_dp, 8.0e-08_dp, failures)
  call check_close('dgenbeta2', dgenbeta2(0.23546097803593721_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), &
    & 1.3948951190910883_dp, 2.0e-08_dp, failures)
  call check_close('pgenbeta2', pgenbeta2(0.23546097803593721_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), &
    & 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varinvbeta', varinvbeta(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.21406440533269799_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dinvbeta', dinvbeta(0.21406440533269799_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 1.0531568747564097_dp, 2.0e-08_dp, failures)
  call check_close('pinvbeta', pinvbeta(0.21406440533269799_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('vargeninvbeta', vargeninvbeta(0.37_dp, a=0.80000000000000004_dp, c=1.2_dp, d=1.5_dp), &
    & 0.39760336259626988_dp, 8.0e-08_dp, failures)
  call check_close('dgeninvbeta', dgeninvbeta(0.39760336259626988_dp, a=0.80000000000000004_dp, c=1.2_dp, &
    & d=1.5_dp), 0.54859056769046721_dp, 2.0e-08_dp, failures)
  call check_close('pgeninvbeta', pgeninvbeta(0.39760336259626988_dp, a=0.80000000000000004_dp, c=1.2_dp, &
    & d=1.5_dp), 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('vartsp', vartsp(0.37_dp, a=1.3_dp, theta=0.40000000000000002_dp), 0.37671695786926401_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dtsp', dtsp(0.37671695786926401_dp, a=1.3_dp, theta=0.40000000000000002_dp), &
    & 1.2768206738570191_dp, 2.0e-08_dp, failures)
  call check_close('ptsp', ptsp(0.37671695786926401_dp, a=1.3_dp, theta=0.40000000000000002_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varkum', varkum(0.37_dp, a=1.2_dp, b=1.5_dp), 0.33075842243724812_dp, 8.0e-08_dp, failures)
  call check_close('dkum', dkum(0.33075842243724812_dp, a=1.2_dp, b=1.5_dp), 1.2367678978575993_dp, 2.0e-08_dp, &
    & failures)
  call check_close('pkum', pkum(0.33075842243724812_dp, a=1.2_dp, b=1.5_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varnormal', varnormal(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & -0.23140935036786164_dp, 8.0e-08_dp, failures)
  call check_close('dnormal', dnormal(-0.23140935036786164_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.29043771861068052_dp, 2.0e-08_dp, failures)
  call check_close('pnormal', pnormal(-0.23140935036786164_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varkumnormal', varkumnormal(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp, b=1.7_dp), -1.059978761056166_dp, 8.0e-08_dp, failures)
  call check_close('dkumnormal', dkumnormal(-1.059978761056166_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.3088601754213266_dp, 2.0e-08_dp, failures)
  call check_close('pkumnormal', pkumnormal(-1.059978761056166_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varexppower', varexppower(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp), -0.1707360819114922_dp, 8.0e-08_dp, failures)
  call check_close('dexppower', dexppower(-0.1707360819114922_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp), 0.28376659718886182_dp, 2.0e-08_dp, failures)
  call check_close('pexppower', pexppower(-0.1707360819114922_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varaep', varaep(0.37_dp, q1=1.3999999999999999_dp, q2=1.8_dp, alpha=0.65000000000000002_dp), &
    & -0.75529325545297143_dp, 8.0e-08_dp, failures)
  call check_close('daep', daep(-0.75529325545297143_dp, q1=1.3999999999999999_dp, q2=1.8_dp, &
    & alpha=0.65000000000000002_dp), 0.3054067002524079_dp, 2.0e-08_dp, failures)
  call check_close('paep', paep(-0.75529325545297143_dp, q1=1.3999999999999999_dp, q2=1.8_dp, &
    & alpha=0.65000000000000002_dp), 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varbetanorm', varbetanorm(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp, b=1.7_dp), -1.0083226131423946_dp, 8.0e-08_dp, failures)
  call check_close('dbetanorm', dbetanorm(-1.0083226131423946_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.30927369383278697_dp, 2.0e-08_dp, failures)
  call check_close('pbetanorm', pbetanorm(-1.0083226131423946_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & a=0.80000000000000004_dp, b=1.7_dp), 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varhalfnorm', varhalfnorm(0.37_dp, sigma=1.3_dp), 0.62624490446014958_dp, 8.0e-08_dp, failures)
  call check_close('dhalfnorm', dhalfnorm(0.62624490446014958_dp, sigma=1.3_dp), 0.54651911441654888_dp, &
    & 2.0e-08_dp, failures)
  call check_close('phalfnorm', phalfnorm(0.62624490446014958_dp, sigma=1.3_dp), 0.37000000000000011_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varkumhalfnorm', varkumhalfnorm(0.37_dp, sigma=1.3_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.27281141827268623_dp, 8.0e-08_dp, failures)
  call check_close('dkumhalfnorm', dkumhalfnorm(0.27281141827268623_dp, sigma=1.3_dp, a=0.80000000000000004_dp, &
    & b=1.7_dp), 0.96652252893656865_dp, 2.0e-08_dp, failures)
  call check_close('pkumhalfnorm', pkumhalfnorm(0.27281141827268623_dp, sigma=1.3_dp, a=0.80000000000000004_dp, &
    & b=1.7_dp), 0.36999999999999977_dp, 2.0e-08_dp, failures)
  call check_close('varT', vart(0.37_dp, n=1.0_dp), -0.43273864224742586_dp, 8.0e-08_dp, failures)
  call check_close('dT', dt(-0.43273864224742586_dp, n=1.0_dp), 0.26810399877969754_dp, 2.0e-08_dp, failures)
  call check_close('pT', pt(-0.43273864224742586_dp, n=1.0_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varast', varast(0.37_dp, nu1=4.5_dp, nu2=6.0_dp, alpha=0.65000000000000002_dp), &
    & -0.7936996520149584_dp, 8.0e-08_dp, failures)
  call check_close('dast', dast(-0.7936996520149584_dp, nu1=4.5_dp, nu2=6.0_dp, alpha=0.65000000000000002_dp), &
    & 0.30417767280930597_dp, 2.0e-08_dp, failures)
  call check_close('past', past(-0.7936996520149584_dp, nu1=4.5_dp, nu2=6.0_dp, alpha=0.65000000000000002_dp), &
    & 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varhalfT', varhalft(0.37_dp, n=1.0_dp), 0.65687722240127944_dp, 8.0e-08_dp, failures)
  call check_close('dhalfT', dhalft(0.65687722240127944_dp, n=1.0_dp), 0.44472598604988028_dp, 2.0e-08_dp, failures)
  call check_close('phalfT', phalft(0.65687722240127944_dp, n=1.0_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varCauchy', varcauchy(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & -0.36256023492165373_dp, 8.0e-08_dp, failures)
  call check_close('dCauchy', dcauchy(-0.36256023492165373_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.2062338452151519_dp, 2.0e-08_dp, failures)
  call check_close('pCauchy', pcauchy(-0.36256023492165373_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varlogcauchy', varlogcauchy(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 24.633880464490872_dp, 8.0e-08_dp, failures)
  call check_close('dlogcauchy', dlogcauchy(24.633880464490872_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.0015677559791456195_dp, 2.0e-08_dp, failures)
  call check_close('plogcauchy', plogcauchy(24.633880464490872_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varhalfcauchy', varhalfcauchy(0.37_dp, sigma=1.3_dp), 0.85394038912166326_dp, 8.0e-08_dp, &
    & failures)
  call check_close('dhalfcauchy', dhalfcauchy(0.85394038912166326_dp, sigma=1.3_dp), 0.34209691234606171_dp, &
    & 2.0e-08_dp, failures)
  call check_close('phalfcauchy', phalfcauchy(0.85394038912166326_dp, sigma=1.3_dp), 0.37000000000000005_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varPCTAlaplace', varpctalaplace(0.37_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp), -3.2555436101482869_dp, 8.0e-08_dp, failures)
  call check_close('dPCTAlaplace', dpctalaplace(-3.2555436101482869_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp), 0.073999999999999982_dp, 2.0e-08_dp, failures)
  call check_close('pPCTAlaplace', ppctalaplace(-3.2555436101482869_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varHBlaplace', varhblaplace(0.37_dp, a=0.80000000000000004_dp, theta=0.59999999999999998_dp, &
    & phi=1.2_dp), -0.042590601691381069_dp, 8.0e-08_dp, failures)
  call check_close('dHBlaplace', dhblaplace(-0.042590601691381069_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp, phi=1.2_dp), 0.44399999999999995_dp, 2.0e-08_dp, failures)
  call check_close('pHBlaplace', phblaplace(-0.042590601691381069_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp, phi=1.2_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varMlaplace', varmlaplace(0.37_dp, theta=0.59999999999999998_dp, phi=1.2_dp, &
    & psi=0.90000000000000002_dp), 0.32900541649447052_dp, 8.0e-08_dp, failures)
  call check_close('dMlaplace', dmlaplace(0.32900541649447052_dp, theta=0.59999999999999998_dp, phi=1.2_dp, &
    & psi=0.90000000000000002_dp), 0.41111111111111115_dp, 2.0e-08_dp, failures)
  call check_close('pMlaplace', pmlaplace(0.32900541649447052_dp, theta=0.59999999999999998_dp, phi=1.2_dp, &
    & psi=0.90000000000000002_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varasylaplace', varasylaplace(0.37_dp, tau=1.2_dp, kappa=1.1000000000000001_dp, &
    & theta=0.59999999999999998_dp), 0.23422656325671654_dp, 8.0e-08_dp, failures)
  call check_close('dasylaplace', dasylaplace(0.23422656325671654_dp, tau=1.2_dp, kappa=1.1000000000000001_dp, &
    & theta=0.59999999999999998_dp), 0.39640834702882222_dp, 2.0e-08_dp, failures)
  call check_close('pasylaplace', pasylaplace(0.23422656325671654_dp, tau=1.2_dp, kappa=1.1000000000000001_dp, &
    & theta=0.59999999999999998_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varlogistic', varlogistic(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & -0.49188185787150068_dp, 8.0e-08_dp, failures)
  call check_close('dlogistic', dlogistic(-0.49188185787150068_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.17930769230769225_dp, 2.0e-08_dp, failures)
  call check_close('plogistic', plogistic(-0.49188185787150068_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varsecant', varsecant(0.37_dp), -0.26754465037494868_dp, 8.0e-08_dp, failures)
  call check_close('dsecant', dsecant(-0.26754465037494868_dp), 0.45887731284199057_dp, 2.0e-08_dp, failures)
  call check_close('psecant', psecant(-0.26754465037494868_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('vargenlogis', vargenlogis(0.37_dp, a=0.80000000000000004_dp, mu=0.20000000000000001_dp, &
    & sigma=1.3_dp), -0.97303703462071822_dp, 8.0e-08_dp, failures)
  call check_close('dgenlogis', dgenlogis(-0.97303703462071822_dp, a=0.80000000000000004_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.16198699085670284_dp, 2.0e-08_dp, failures)
  call check_close('pgenlogis', pgenlogis(-0.97303703462071822_dp, a=0.80000000000000004_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varhalflogis', varhalflogis(0.37_dp, lambda=1.3_dp), 0.59757399956660939_dp, 8.0e-08_dp, failures)
  call check_close('dhalflogis', dhalflogis(0.59757399956660939_dp, lambda=1.3_dp), 0.56101499999999993_dp, &
    & 2.0e-08_dp, failures)
  call check_close('phalflogis', phalflogis(0.59757399956660939_dp, lambda=1.3_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varloglogis', varloglogis(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.58495973376908295_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dloglogis', dloglogis(0.58495973376908295_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.67743124376562736_dp, 2.0e-08_dp, failures)
  call check_close('ploglogis', ploglogis(0.58495973376908295_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varexplogis', varexplogis(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), -1.5339715068117081_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dexplogis', dexplogis(-1.5339715068117081_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.12387240477277278_dp, 2.0e-08_dp, failures)
  call check_close('pexplogis', pexplogis(-1.5339715068117081_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varHlogis', varhlogis(0.37_dp, k=0.40000000000000002_dp), -0.59311123763200579_dp, 8.0e-08_dp, &
    & failures)
  call check_close('dHlogis', dhlogis(-0.59311123763200579_dp, k=0.40000000000000002_dp), 0.18840253557972139_dp, &
    & 2.0e-08_dp, failures)
  call check_close('pHlogis', phlogis(-0.59311123763200579_dp, k=0.40000000000000002_dp), 0.37_dp, 2.0e-08_dp, &
    & failures)
  call check_close('varlognorm', varlognorm(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.79341461498631927_dp, 8.0e-08_dp, failures)
  call check_close('dlognorm', dlognorm(0.79341461498631927_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.36606045959424188_dp, 2.0e-08_dp, failures)
  call check_close('plognorm', plognorm(0.79341461498631927_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varbetalognorm', varbetalognorm(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.36483042837587559_dp, 8.0e-08_dp, failures)
  call check_close('dbetalognorm', dbetalognorm(0.36483042837587559_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.84771902170986158_dp, 2.0e-08_dp, failures)
  call check_close('pbetalognorm', pbetalognorm(0.36483042837587559_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varburr', varburr(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.58495973376908295_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dburr', dburr(0.58495973376908295_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.67743124376562724_dp, 2.0e-08_dp, failures)
  call check_close('pburr', pburr(0.58495973376908295_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, 2.0e-08_dp, &
    & failures)
  call check_close('varbetaburr', varbetaburr(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp), &
    & 0.51831752389781294_dp, 8.0e-08_dp, failures)
  call check_close('dbetaburr', dbetaburr(0.51831752389781294_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), 0.89425510277812803_dp, 2.0e-08_dp, failures)
  call check_close('pbetaburr', pbetaburr(0.51831752389781294_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varkumburr7', varkumburr7(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & k=0.40000000000000002_dp, c=1.2_dp), 0.63085227808526534_dp, 8.0e-08_dp, failures)
  call check_close('dkumburr7', dkumburr7(0.63085227808526534_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & k=0.40000000000000002_dp, c=1.2_dp), 0.37298107777954553_dp, 2.0e-08_dp, failures)
  call check_close('pkumburr7', pkumburr7(0.63085227808526534_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & k=0.40000000000000002_dp, c=1.2_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varbetaburr7', varbetaburr7(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & k=0.40000000000000002_dp), 0.67508764486614781_dp, 8.0e-08_dp, failures)
  call check_close('dbetaburr7', dbetaburr7(0.67508764486614781_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & k=0.40000000000000002_dp), 0.34933632470985743_dp, 2.0e-08_dp, failures)
  call check_close('pbetaburr7', pbetaburr7(0.67508764486614781_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & k=0.40000000000000002_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varlomax', varlomax(0.37_dp, a=0.80000000000000004_dp, lambda=1.3_dp), 1.0161552449515709_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dlomax', dlomax(1.0161552449515709_dp, a=0.80000000000000004_dp, lambda=1.3_dp), &
    & 0.21760199412303999_dp, 2.0e-08_dp, failures)
  call check_close('plomax', plomax(1.0161552449515709_dp, a=0.80000000000000004_dp, lambda=1.3_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varbetalomax', varbetalomax(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, lambda=1.3_dp), 0.45204803286602047_dp, 8.0e-08_dp, failures)
  call check_close('dbetalomax', dbetalomax(0.45204803286602047_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, lambda=1.3_dp), 0.47435353544963182_dp, 2.0e-08_dp, failures)
  call check_close('pbetalomax', pbetalomax(0.45204803286602047_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, lambda=1.3_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('vargumbel', vargumbel(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.20749360092748778_dp, 8.0e-08_dp, failures)
  call check_close('dgumbel', dgumbel(0.20749360092748778_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & 0.28297949318248516_dp, 2.0e-08_dp, failures)
  call check_close('pgumbel', pgumbel(0.20749360092748778_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varkumgumbel', varkumgumbel(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), -0.56010701234295079_dp, 8.0e-08_dp, failures)
  call check_close('dkumgumbel', dkumgumbel(-0.56010701234295079_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.36935740955858354_dp, 2.0e-08_dp, failures)
  call check_close('pkumgumbel', pkumgumbel(-0.56010701234295079_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varbetagumbel', varbetagumbel(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), -0.51664833743164662_dp, 8.0e-08_dp, failures)
  call check_close('dbetagumbel', dbetagumbel(-0.51664833743164662_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.36538316167096074_dp, 2.0e-08_dp, failures)
  call check_close('pbetagumbel', pbetagumbel(-0.51664833743164662_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & mu=0.20000000000000001_dp, sigma=1.3_dp), 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('vargumbel2', vargumbel2(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 5.0958544789085529_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dgumbel2', dgumbel2(5.0958544789085529_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.045697119609769861_dp, 2.0e-08_dp, failures)
  call check_close('pgumbel2', pgumbel2(5.0958544789085529_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.37_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varbetagumbel2', varbetagumbel2(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp, &
    & d=1.5_dp), 6.2825260698102747_dp, 8.0e-08_dp, failures)
  call check_close('dbetagumbel2', dbetagumbel2(6.2825260698102747_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & c=1.2_dp, d=1.5_dp), 0.041943310259904353_dp, 2.0e-08_dp, failures)
  call check_close('pbetagumbel2', pbetagumbel2(6.2825260698102747_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & c=1.2_dp, d=1.5_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varfrechet', varfrechet(0.37_dp, alpha=0.65000000000000002_dp, sigma=1.3_dp), &
    & 1.3115798871055921_dp, 8.0e-08_dp, failures)
  call check_close('dfrechet', dfrechet(1.3115798871055921_dp, alpha=0.65000000000000002_dp, sigma=1.3_dp), &
    & 0.18231270095707802_dp, 2.0e-08_dp, failures)
  call check_close('pfrechet', pfrechet(1.3115798871055921_dp, alpha=0.65000000000000002_dp, sigma=1.3_dp), &
    & 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varbetafrechet', varbetafrechet(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.55669360183965433_dp, 8.0e-08_dp, failures)
  call check_close('dbetafrechet', dbetafrechet(0.55669360183965433_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.55461167613866613_dp, 2.0e-08_dp, failures)
  call check_close('pbetafrechet', pbetafrechet(0.55669360183965433_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.37000000000000038_dp, 2.0e-08_dp, failures)
  call check_close('varWeibull', varweibull(0.37_dp, alpha=0.65000000000000002_dp, sigma=1.3_dp), &
    & 0.39633206479755245_dp, 8.0e-08_dp, failures)
  call check_close('dWeibull', dweibull(0.39633206479755245_dp, alpha=0.65000000000000002_dp, sigma=1.3_dp), &
    & 0.47738635732497831_dp, 2.0e-08_dp, failures)
  call check_close('pWeibull', pweibull(0.39633206479755245_dp, alpha=0.65000000000000002_dp, sigma=1.3_dp), &
    & 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varlogisrayleigh', varlogisrayleigh(0.37_dp, a=0.80000000000000004_dp, lambda=1.3_dp), &
    & 0.79888746718284109_dp, 8.0e-08_dp, failures)
  call check_close('dlogisrayleigh', dlogisrayleigh(0.79888746718284109_dp, a=0.80000000000000004_dp, &
    & lambda=1.3_dp), 0.57036006909929471_dp, 2.0e-08_dp, failures)
  call check_close('plogisrayleigh', plogisrayleigh(0.79888746718284109_dp, a=0.80000000000000004_dp, &
    & lambda=1.3_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varmoweibull', varmoweibull(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, lambda=1.3_dp), &
    & 0.55873754525171981_dp, 8.0e-08_dp, failures)
  call check_close('dmoweibull', dmoweibull(0.55873754525171981_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & lambda=1.3_dp), 1.2533635458842014_dp, 2.0e-08_dp, failures)
  call check_close('pmoweibull', pmoweibull(0.55873754525171981_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & lambda=1.3_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varbetaweibull', varbetaweibull(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.10427093572912915_dp, 8.0e-08_dp, failures)
  call check_close('dbetaweibull', dbetaweibull(0.10427093572912915_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 1.5460652762993374_dp, 2.0e-08_dp, failures)
  call check_close('pbetaweibull', pbetaweibull(0.10427093572912915_dp, a=0.80000000000000004_dp, b=1.7_dp, &
    & alpha=0.65000000000000002_dp, sigma=1.3_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('vardweibull', vardweibull(0.37_dp, c=1.2_dp, mu=0.20000000000000001_dp, sigma=1.3_dp), &
    & -0.27812534940342765_dp, 8.0e-08_dp, failures)
  call check_close('ddweibull', ddweibull(-0.27812534940342765_dp, c=1.2_dp, mu=0.20000000000000001_dp, &
    & sigma=1.3_dp), 0.0_dp, 2.0e-08_dp, failures)
  call check_close('pdweibull', pdweibull(-0.27812534940342765_dp, c=1.2_dp, mu=0.20000000000000001_dp, &
    & sigma=1.3_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('vargenpowerweibull', vargenpowerweibull(0.37_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp), 0.85637274487360426_dp, 8.0e-08_dp, failures)
  call check_close('dgenpowerweibull', dgenpowerweibull(0.85637274487360426_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp), 0.24214599041364138_dp, 2.0e-08_dp, failures)
  call check_close('pgenpowerweibull', pgenpowerweibull(0.85637274487360426_dp, a=0.80000000000000004_dp, &
    & theta=0.59999999999999998_dp), 0.37000000000000011_dp, 2.0e-08_dp, failures)
  call check_close('varchen', varchen(0.37_dp, b=1.7_dp, lambda=1.3_dp), 0.49647403198040407_dp, 8.0e-08_dp, failures)
  call check_close('dchen', dchen(0.49647403198040407_dp, b=1.7_dp, lambda=1.3_dp), 1.1559304153366319_dp, &
    & 2.0e-08_dp, failures)
  call check_close('pchen', pchen(0.49647403198040407_dp, b=1.7_dp, lambda=1.3_dp), 0.36999999999999988_dp, &
    & 2.0e-08_dp, failures)
  call check_close('varxie', varxie(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, lambda=1.3_dp), &
    & 0.44404688667105785_dp, 8.0e-08_dp, failures)
  call check_close('dxie', dxie(0.44404688667105785_dp, a=0.80000000000000004_dp, b=1.7_dp, lambda=1.3_dp), &
    & 1.331732072290075_dp, 2.0e-08_dp, failures)
  call check_close('pxie', pxie(0.44404688667105785_dp, a=0.80000000000000004_dp, b=1.7_dp, lambda=1.3_dp), &
    & 0.36999999999999988_dp, 2.0e-08_dp, failures)
  call check_close('varTL', vartl(0.37_dp, lambda=1.3_dp), -0.21067798087476547_dp, 8.0e-08_dp, failures)
  call check_close('varRS', varrs(0.37_dp, b=1.7_dp, c=1.2_dp, d=1.5_dp), -0.25994337819899432_dp, 8.0e-08_dp, &
    & failures)
  call check_close('varFR', varfr(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), -0.15630775969280442_dp, &
    & 8.0e-08_dp, failures)
  call check_close('varHL', varhl(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp, c=1.2_dp), 1.1881292596406579_dp, &
    & 8.0e-08_dp, failures)
  call check_close('varloglog', varloglog(0.37_dp, a=0.80000000000000004_dp, lambda=1.3_dp), 1.5880171114101744_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dloglog', dloglog(1.5880171114101744_dp, a=0.80000000000000004_dp, lambda=1.3_dp), &
    & 0.1762471439818343_dp, 2.0e-08_dp, failures)
  call check_close('ploglog', ploglog(1.5880171114101744_dp, a=0.80000000000000004_dp, lambda=1.3_dp), &
    & 0.37000000000000011_dp, 2.0e-08_dp, failures)
  call check_close('varexplog', varexplog(0.37_dp, a=0.80000000000000004_dp, b=1.7_dp), 0.2482379033440672_dp, &
    & 8.0e-08_dp, failures)
  call check_close('dexplog', dexplog(0.2482379033440672_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 1.1499360441322648_dp, 2.0e-08_dp, failures)
  call check_close('pexplog', pexplog(0.2482379033440672_dp, a=0.80000000000000004_dp, b=1.7_dp), &
    & 0.37000000000000011_dp, 2.0e-08_dp, failures)
  call check_close('varexpgeo', varexpgeo(0.37_dp, theta=0.59999999999999998_dp, lambda=1.3_dp), &
    & 0.54146151718169777_dp, 8.0e-08_dp, failures)
  call check_close('dexpgeo', dexpgeo(0.54146151718169777_dp, theta=0.59999999999999998_dp, lambda=1.3_dp), &
    & 0.59964666666666666_dp, 2.0e-08_dp, failures)
  call check_close('pexpgeo', pexpgeo(0.54146151718169777_dp, theta=0.59999999999999998_dp, lambda=1.3_dp), &
    & 0.36999999999999994_dp, 2.0e-08_dp, failures)
  call check_close('varexppois', varexppois(0.37_dp, b=1.7_dp, lambda=1.3_dp), 0.16236603632319493_dp, 8.0e-08_dp, &
    & failures)
  call check_close('dexppois', dexppois(0.16236603632319493_dp, b=1.7_dp, lambda=1.3_dp), 1.6847043465422777_dp, &
    & 2.0e-08_dp, failures)
  call check_close('pexppois', pexppois(0.16236603632319493_dp, b=1.7_dp, lambda=1.3_dp), 0.37_dp, 2.0e-08_dp, &
    & failures)
  call check_close('varTL2', vartl2(0.37_dp, b=1.7_dp), 0.33455807921659497_dp, 8.0e-08_dp, failures)
  call check_close('dTL2', dtl2(0.33455807921659497_dp, b=1.7_dp), 1.5024145594360636_dp, 2.0e-08_dp, failures)
  call check_close('pTL2', ptl2(0.33455807921659497_dp, b=1.7_dp), 0.37_dp, 2.0e-08_dp, failures)
  call check_close('varquad', varquad(0.37_dp, a=0.0_dp, b=2.0_dp), 0.3617495701140091_dp, 8.0e-08_dp, failures)
  call check_close('dquad', dquad(0.3617495701140091_dp, a=0.0_dp, b=2.0_dp), 0.61104541687447833_dp, 2.0e-08_dp, &
    & failures)
  call check_close('pquad', pquad(0.3617495701140091_dp, a=0.0_dp, b=2.0_dp), 0.36999999999999988_dp, 2.0e-08_dp, &
    & failures)
  call check_close('varschabe', varschabe(0.37_dp, gamma=0.69999999999999996_dp, theta=0.59999999999999998_dp), &
    & 0.1168421052631579_dp, 8.0e-08_dp, failures)
  call check_close('dschabe', dschabe(0.1168421052631579_dp, gamma=0.69999999999999996_dp, &
    & theta=0.59999999999999998_dp), 3.0362745098039214_dp, 2.0e-08_dp, failures)
  call check_close('pschabe', pschabe(0.1168421052631579_dp, gamma=0.69999999999999996_dp, &
    & theta=0.59999999999999998_dp), 0.37000000000000005_dp, 2.0e-08_dp, failures)
  call check_close('varBS', varbs(0.37_dp, gamma=0.69999999999999996_dp), 0.79312202914662233_dp, 8.0e-08_dp, &
    & failures)
  call check_close('dBS', dbs(0.79312202914662233_dp, gamma=0.69999999999999996_dp), 0.68464928726143337_dp, &
    & 2.0e-08_dp, failures)
  call check_close('pBS', pbs(0.79312202914662233_dp, gamma=0.69999999999999996_dp), 0.36999999999999977_dp, &
    & 2.0e-08_dp, failures)
  call check_close('vargev', vargev(0.37_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, xi=0.20000000000000001_dp), &
    & 0.20749792213059237_dp, 8.0e-08_dp, failures)
  call check_close('dgev', dgev(0.20749792213059237_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & xi=0.20000000000000001_dp), 0.28265344494861311_dp, 2.0e-08_dp, failures)
  call check_close('pgev', pgev(0.20749792213059237_dp, mu=0.20000000000000001_dp, sigma=1.3_dp, &
    & xi=0.20000000000000001_dp), 0.36999999999999983_dp, 2.0e-08_dp, failures)
  if (failures /= 0) error stop 'VaRES reference tests failed'
  print '(a,i0,a)', 'VaRES source-reference checks passed: ', 292, '.'
contains
  subroutine check_close(label,actual,expected,tol,failures)
    character(*),intent(in)::label
    real(dp),intent(in)::actual,expected,tol
    integer,intent(inout)::failures
    if(.not.ieee_is_finite(actual).or.abs(actual-expected)>tol*(1.0_dp+abs(expected)))then
      failures=failures+1
      print '(a,2es24.14)',trim(label)//' failed: ',actual,expected
    end if
  end subroutine check_close
end program test_reference
