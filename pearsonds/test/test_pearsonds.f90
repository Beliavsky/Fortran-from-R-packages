program test_pearsonds
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf
use pearsonds_mod
implicit none
integer :: failures, fam
real(kind=dp) :: tol, x, p, q, m(4), m2(4)
type(pearson_params_t) :: fit
type(pearson_ml_result_t) :: ml
complex(kind=dp) :: hg
real(kind=dp), allocatable :: v(:)

failures = 0
tol = 2.0e-6_dp

call check(abs(dpearson0(0.0_dp,0.0_dp,1.0_dp)-0.3989422804014327_dp)<1.0e-12_dp,"normal density")
call check(abs(ppearson0(0.0_dp,0.0_dp,1.0_dp)-0.5_dp)<1.0e-12_dp,"normal cdf")
call check(abs(qpearson0(0.5_dp,0.0_dp,1.0_dp))<1.0e-10_dp,"normal quantile")

x = 1.1_dp
p = ppearsoni(x,2.5_dp,4.0_dp,-1.0_dp,5.0_dp)
q = qpearsoni(p,2.5_dp,4.0_dp,-1.0_dp,5.0_dp)
call check(abs(q-x)<tol,"type I p/q inverse")
p = ppearsonii(x,3.0_dp,-2.0_dp,6.0_dp)
q = qpearsonii(p,3.0_dp,-2.0_dp,6.0_dp)
call check(abs(q-x)<tol,"type II p/q inverse")
p = ppearsoniii(x,4.0_dp,-1.0_dp,2.0_dp)
q = qpearsoniii(p,4.0_dp,-1.0_dp,2.0_dp)
call check(abs(q-x)<tol,"type III p/q inverse")
p = ppearsonv(3.0_dp,6.0_dp,0.0_dp,2.0_dp)
q = qpearsonv(p,6.0_dp,0.0_dp,2.0_dp)
call check(abs(q-3.0_dp)<tol,"type V p/q inverse")
p = ppearsonvi(3.0_dp,4.0_dp,7.0_dp,0.0_dp,2.0_dp)
q = qpearsonvi(p,4.0_dp,7.0_dp,0.0_dp,2.0_dp)
call check(abs(q-3.0_dp)<2.0e-5_dp,"type VI p/q inverse")
p = ppearsonvii(1.0_dp,8.0_dp,0.0_dp,1.5_dp)
q = qpearsonvii(p,8.0_dp,0.0_dp,1.5_dp)
call check(abs(q-1.0_dp)<2.0e-5_dp,"type VII p/q inverse")

! Type IV: CDF is the same numerical-integration fallback used by PearsonDS when GSL is absent.
p = ppearsoniv(0.25_dp,5.0_dp,1.25_dp,-0.2_dp,1.4_dp,tol=1.0e-7_dp)
q = qpearsoniv(p,5.0_dp,1.25_dp,-0.2_dp,1.4_dp,tol=2.0e-6_dp)
call check(p>0.0_dp.and.p<1.0_dp.and.abs(q-0.25_dp)<2.0e-3_dp,"type IV p/q inverse")
call check(ieee_is_finite(log_pearson_iv_norm(5.0_dp,1.25_dp,1.4_dp)),"type IV log normalization")
call check(ppearsoniv(-ieee_value(0.0_dp,ieee_positive_inf),5.0_dp,1.25_dp,-0.2_dp,1.4_dp)==0.0_dp,"type IV negative infinity")
call check(ppearsoniv( ieee_value(0.0_dp,ieee_positive_inf),5.0_dp,1.25_dp,-0.2_dp,1.4_dp)==1.0_dp,"type IV positive infinity")
hg=hypergeom_2f1(cmplx(1.0_dp,0.0_dp,dp),cmplx(1.0_dp,0.0_dp,dp), &
                 cmplx(2.0_dp,0.0_dp,dp),cmplx(0.2_dp,0.0_dp,dp))
call check(abs(real(hg,dp)+log(0.8_dp)/0.2_dp)<1.0e-10_dp,"2F1 series")

m = pearson0moments(1.0_dp,2.0_dp)
fit = pearson_fit_m(m(1),m(2),m(3),m(4))
call check(fit%family==0,"moment classifier type 0")

m = pearsonimoments(2.0_dp,5.0_dp,-1.0_dp,4.0_dp)
fit = pearson_fit_m(m(1),m(2),m(3),m(4))
call check(fit%family==1,"moment classifier type I")
if(fit%family==1)then
 m2=pearson_moments(fit)
 call check(maxval(abs(m2-m))<2.0e-6_dp,"type I moment roundtrip")
end if

m = pearsoniimoments(3.0_dp,-2.0_dp,6.0_dp)
fit = pearson_fit_m(m(1),m(2),m(3),m(4))
call check(fit%family==2,"moment classifier type II")

m = pearsoniiimoments(5.0_dp,-1.0_dp,2.0_dp)
fit = pearson_fit_m(m(1),m(2),m(3),m(4))
call check(fit%family==3,"moment classifier type III")

m = pearsonivmoments(6.0_dp,2.0_dp,0.5_dp,1.25_dp)
fit = pearson_fit_m(m(1),m(2),m(3),m(4))
call check(fit%family==4,"moment classifier type IV")
if(fit%family==4)then
 m2=pearson_moments(fit)
 call check(maxval(abs(m2-m))<5.0e-6_dp,"type IV moment roundtrip")
end if

m = pearsonvmoments(8.0_dp,0.0_dp,2.0_dp)
fit = pearson_fit_m(m(1),m(2),m(3),m(4))
call check(fit%family==5,"moment classifier type V")

m = pearsonvimoments(3.0_dp,8.0_dp,-1.0_dp,2.0_dp)
fit = pearson_fit_m(m(1),m(2),m(3),m(4))
call check(fit%family==6,"moment classifier type VI")

m = pearsonviimoments(9.0_dp,0.5_dp,1.25_dp)
fit = pearson_fit_m(m(1),m(2),m(3),m(4))
call check(fit%family==7,"moment classifier type VII")

! Maximum-likelihood exact normal fit.
ml=pearson_fit_ml_type([-1.0_dp,0.0_dp,1.0_dp,2.0_dp],0)
call check(ml%convergence==0.and.abs(ml%fit%par(1)-0.5_dp)<1.0e-12_dp,"normal MLE")

! Generic dispatch smoke test.
fit%family=3
fit%npar=3
fit%par=0.0_dp
fit%par(1:3)=[5.0_dp,-1.0_dp,2.0_dp]
call check(dpearson(0.0_dp,fit)>0.0_dp,"generic scalar density dispatch")
v=dpearson([0.0_dp,1.0_dp,2.0_dp],fit)
call check(size(v)==3.and.all(v>=0.0_dp),"generic density dispatch")

! RNG smoke tests for all families except IV m=1 edge case.
do fam=0,7
 select case(fam)
 case(0)
 fit%family=0
 fit%npar=2
 fit%par(1:2)=[0.0_dp,1.0_dp]
 case(1)
 fit%family=1
 fit%npar=4
 fit%par=[2.0_dp,3.0_dp,0.0_dp,1.0_dp]
 case(2)
 fit%family=2
 fit%npar=3
 fit%par(1:3)=[2.0_dp,0.0_dp,1.0_dp]
 case(3)
 fit%family=3
 fit%npar=3
 fit%par(1:3)=[3.0_dp,0.0_dp,1.0_dp]
 case(4)
 fit%family=4
 fit%npar=4
 fit%par=[4.0_dp,1.0_dp,0.0_dp,1.0_dp]
 case(5)
 fit%family=5
 fit%npar=3
 fit%par(1:3)=[6.0_dp,0.0_dp,1.0_dp]
 case(6)
 fit%family=6
 fit%npar=4
 fit%par=[3.0_dp,7.0_dp,0.0_dp,1.0_dp]
 case(7)
 fit%family=7
 fit%npar=3
 fit%par(1:3)=[8.0_dp,0.0_dp,1.0_dp]
 end select
 v=rpearson(5,fit)
 call check(size(v)==5.and.all(ieee_is_finite(v)),"rng family")
end do

if(failures>0)then
 write(*,'(a,i0)') 'FAILURES: ',failures
 error stop 1
else
 write(*,'(a)') 'All PearsonDS translation tests passed.'
end if

contains
subroutine check(ok,name)
logical,intent(in)::ok
character(len=*),intent(in)::name
if(.not.ok)then
 failures=failures+1
 write(*,'(a)') 'FAIL: '//trim(name)
end if
end subroutine check
end program test_pearsonds
