! SPDX-License-Identifier: GPL-3.0-only
program test_fft_stress
  use poisson_binomial, only : dp, dpb_convolve, dpb_dividefft
  implicit none
  real(dp)::p(96)
  real(dp),allocatable::a(:),b(:)
  integer::i
  do i=1,size(p)
    p(i)=0.02_dp+0.96_dp*real(mod(37*i,97),dp)/96.0_dp
  end do
  a=dpb_convolve(p); b=dpb_dividefft(p)
  call chk(maxval(abs(a-b))<2e-12_dp,"FFT/convolution agreement")
  call chk(abs(sum(b)-1.0_dp)<2e-14_dp,"FFT normalization")
  print '(a)', 'test_fft_stress: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_fft_stress
