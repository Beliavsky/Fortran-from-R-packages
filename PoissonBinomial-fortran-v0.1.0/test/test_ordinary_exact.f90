! SPDX-License-Identifier: GPL-3.0-only
program test_ordinary_exact
  use poisson_binomial, only : dp, dpb_convolve, dpb_dividefft, dpb_characteristic, &
                               dpb_recursive
  implicit none
  real(dp), parameter :: probs(4) = [0.1_dp,0.3_dp,0.6_dp,0.8_dp]
  real(dp), parameter :: ref(0:4) = [0.0504_dp,0.3044_dp,0.4544_dp,0.1764_dp,0.0144_dp]
  real(dp), allocatable :: a(:),b(:),c(:),d(:)
  a=dpb_convolve(probs)
  b=dpb_dividefft(probs)
  c=dpb_characteristic(probs)
  d=dpb_recursive(probs)
  call chk(maxval(abs(a-ref))<2.0e-14_dp,"convolve reference")
  call chk(maxval(abs(b-ref))<2.0e-14_dp,"divide reference")
  call chk(maxval(abs(c-ref))<2.0e-13_dp,"characteristic reference")
  call chk(maxval(abs(d-ref))<2.0e-13_dp,"recursive reference")
  print '(a)', 'test_ordinary_exact: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_ordinary_exact
