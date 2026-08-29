program test_compquadform
use r_compat, only: dp
use compquadform_mod, only: davies_result_t, farebrother_result_t, &
   imhof_result_t, davies, farebrother, imhof, liu
implicit none
integer :: failures
failures = 0
call test_davies(failures)
call test_farebrother(failures)
call test_imhof(failures)
call test_liu(failures)
if (failures /= 0) then
   write(*,'(a,i0)') 'FAIL: ', failures
   error stop 1
end if
write(*,'(a)') 'All CompQuadForm translation tests passed.'

contains

subroutine test_davies(nfail)
integer, intent(inout) :: nfail
type(davies_result_t) :: fit
real(kind=dp) :: lambda(3), delta(3)
integer :: h(3)
lambda = [1.0_dp, 0.5_dp, -0.25_dp]
delta = [0.2_dp, 1.0_dp, 0.0_dp]
h = [2, 3, 1]
fit = davies(2.0_dp, lambda, h, delta, sigma=0.3_dp, &
   lim=100000, acc=1.0e-6_dp)
call assert_int('davies general ifault', fit%ifault, 0, nfail)
call assert_close('davies general Qq', fit%qq, &
   0.74590020254474165_dp, 2.0e-13_dp, nfail)
call assert_close('davies trace absolute sum', fit%trace(1), &
   1.1774947459686436_dp, 3.0e-13_dp, nfail)
call assert_close('davies trace terms', fit%trace(2), 67.0_dp, 0.0_dp, nfail)

fit = davies(2.0_dp, lambda, h, delta, sigma=0.0_dp, &
   lim=100000, acc=1.0e-8_dp)
call assert_int('davies mixed-sign ifault', fit%ifault, 0, nfail)
call assert_close('davies mixed-sign Qq', fit%qq, &
   0.74657588081445292_dp, 3.0e-13_dp, nfail)
end subroutine test_davies

subroutine test_farebrother(nfail)
integer, intent(inout) :: nfail
type(farebrother_result_t) :: fit
real(kind=dp) :: lambda(3), delta(3)
integer :: h(3)
lambda = [0.5_dp, 1.2_dp, 2.0_dp]
delta = [0.0_dp, 0.7_dp, 0.3_dp]
h = [1, 2, 3]
fit = farebrother(6.0_dp, lambda, h, delta, maxit=100000, &
   eps=1.0e-12_dp, mode=1.0_dp)
call assert_int('farebrother general ifault', fit%ifault, 0, nfail)
call assert_close('farebrother general Qq', fit%qq, &
   0.73174782484189593_dp, 3.0e-13_dp, nfail)
call assert_close('farebrother density', fit%dnsty, &
   0.077123912630009214_dp, 3.0e-14_dp, nfail)
end subroutine test_farebrother

subroutine test_imhof(nfail)
integer, intent(inout) :: nfail
type(imhof_result_t) :: fit
real(kind=dp) :: lambda(3), delta(3)
integer :: h(3)
lambda = [1.0_dp, 0.5_dp, -0.25_dp]
delta = [0.2_dp, 1.0_dp, 0.0_dp]
h = [2, 3, 1]
fit = imhof(2.0_dp, lambda, h, delta, epsabs=1.0e-8_dp, &
   epsrel=1.0e-8_dp, limit=10000)
call assert_close('imhof mixed-sign Qq', fit%qq, &
   0.7465758808999261_dp, 2.0e-8_dp, nfail)

lambda = [0.5_dp, 1.2_dp, 2.0_dp]
delta = [0.0_dp, 0.7_dp, 0.3_dp]
h = [1, 2, 3]
fit = imhof(6.0_dp, lambda, h, delta, epsabs=1.0e-8_dp, &
   epsrel=1.0e-8_dp, limit=10000)
call assert_close('imhof positive Qq', fit%qq, &
   0.7317478248427531_dp, 2.0e-8_dp, nfail)
end subroutine test_imhof

subroutine test_liu(nfail)
integer, intent(inout) :: nfail
real(kind=dp) :: lambda(3), delta(3), qq
integer :: h(3)
lambda = [0.5_dp, 1.2_dp, 2.0_dp]
delta = [0.0_dp, 0.7_dp, 0.3_dp]
h = [1, 2, 3]
qq = liu(6.0_dp, lambda, h, delta)
call assert_close('liu general Qq', qq, 0.730099834902864_dp, &
   2.0e-12_dp, nfail)

qq = liu(7.0_dp, [1.0_dp], [5], [1.2_dp])
call assert_close('liu single noncentral chi-square', qq, &
   0.3440931505495764_dp, 2.0e-12_dp, nfail)
end subroutine test_liu

subroutine assert_close(name, actual, expected, tol, nfail)
character(len=*), intent(in) :: name
real(kind=dp), intent(in) :: actual, expected, tol
integer, intent(inout) :: nfail
if (abs(actual - expected) > tol) then
   write(*,'(a,2(1x,es24.16),a,es12.4)') trim(name)//':', &
      actual, expected, ' tol=', tol
   nfail = nfail + 1
end if
end subroutine assert_close

subroutine assert_int(name, actual, expected, nfail)
character(len=*), intent(in) :: name
integer, intent(in) :: actual, expected
integer, intent(inout) :: nfail
if (actual /= expected) then
   write(*,'(a,2(1x,i0))') trim(name)//':', actual, expected
   nfail = nfail + 1
end if
end subroutine assert_int

end program test_compquadform
