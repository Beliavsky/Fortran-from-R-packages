program test_series_fourier
use tweedie
implicit none
real(dp), parameter :: yref(8) = [0.001_dp,0.01_dp,0.05_dp,0.1_dp,0.5_dp,0.75_dp,1.0_dp,2.0_dp]
real(dp), parameter :: dref(8) = [ &
 0.15157581294443723_dp,0.15106506794690938_dp,0.14881507399160360_dp, &
 0.14604794879022837_dp,0.12562774249866016_dp,0.11429584298035371_dp, &
 0.10395520767485422_dp,0.070939964617860465_dp]
real(dp), parameter :: pref(8) = [ &
 0.60668226394985880_dp,0.60804414666879403_dp,0.61404164111086212_dp, &
 0.62141300819991796_dp,0.67564929629490433_dp,0.70561818350196892_dp, &
 0.73287980379682038_dp,0.81930997272516115_dp]
real(dp)::ds,di,ps,pii,dd
integer::i,failures,es,regions
real(dp)::relerr
failures=0

! Independent compound-Poisson/Gamma references for p=1.5, mu=1, phi=4.
do i=1,size(yref)
 ds=dtweedie_series(yref(i),1.5_dp,1.0_dp,4.0_dp)
 ps=ptweedie_series(yref(i),1.5_dp,1.0_dp,4.0_dp)
 call check_close('series density',ds,dref(i),2.0e-12_dp,failures)
 call check_close('series cdf',ps,pref(i),2.0e-12_dp,failures)
end do

! Fourier inversion agrees with the series in both Tweedie regions.
di=dtweedie_inversion(1.0_dp,1.0_dp,4.0_dp,1.5_dp,exitstatus=es,relerr=relerr,regions=regions)
call check_close('small-p inversion density',di,dref(7),2.0e-10_dp,failures)
pii=ptweedie_inversion(1.0_dp,1.0_dp,4.0_dp,1.5_dp,exitstatus=es,relerr=relerr,regions=regions)
call check_close('small-p inversion cdf',pii,pref(7),2.0e-10_dp,failures)

ds=dtweedie_series(1.0_dp,2.5_dp,1.3_dp,0.8_dp)
di=dtweedie_inversion(1.0_dp,1.3_dp,0.8_dp,2.5_dp,exitstatus=es,relerr=relerr,regions=regions)
call check_close('big-p series/inversion',di,ds,2.0e-7_dp,failures)

! The stored interpolation grid should track direct inversion closely.
dd=dtweedie(1.0_dp,1.0_dp,0.4_dp,1.5_dp)
di=dtweedie_inversion(1.0_dp,1.0_dp,0.4_dp,1.5_dp)
call check_close('stored interpolation',dd,di,2.0e-5_dp,failures)

! Force Fourier inversion for p=3 and compare the exact inverse-Gaussian result.
di=dtweedie_inversion(1.0_dp,1.4_dp,0.74_dp,3.0_dp,igexact=.false.)
call check_close('p=3 forced inversion',di,0.4388738851125460_dp,2.0e-10_dp,failures)
pii=ptweedie_inversion(1.0_dp,1.4_dp,0.74_dp,3.0_dp,igexact=.false.)
call check_close('p=3 forced cdf inversion',pii,0.5294018089021599_dp,2.0e-10_dp,failures)

if(failures/=0)error stop 'test_series_fourier failed'
print '(a)','test_series_fourier: PASS'
contains
subroutine check_close(name,got,expected,tol,failures)
character(*),intent(in)::name
real(dp),intent(in)::got,expected,tol
integer,intent(inout)::failures
if(abs(got-expected)>tol*max(1.0_dp,abs(expected)))then
 print '(a,2es24.14)',trim(name)//' got/expected: ',got,expected
 failures=failures+1
end if
end subroutine
end program test_series_fourier
