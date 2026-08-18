! SPDX-License-Identifier: GPL-2.0-or-later
module test_hitro_target
  use runuran_kinds, only : dp, pi
  implicit none
contains
  real(dp) function bivnorm_pdf(x,p) result(y)
    real(dp),intent(in)::x(:),p(:)
    if(size(p)<0) stop
    y=exp(-0.5_dp*sum(x*x))/(2.0_dp*pi)
  end function bivnorm_pdf
  real(dp) function bivnorm_logpdf(x,p) result(y)
    real(dp),intent(in)::x(:),p(:)
    if(size(p)<0) stop
    y=-0.5_dp*sum(x*x)-log(2.0_dp*pi)
  end function bivnorm_logpdf
end module test_hitro_target

program test_hitro
  use runuran
  use test_hitro_target
  implicit none
  type(rng_state)::rng
  type(multivariate_distribution)::d
  type(mv_generator)::g
  real(dp)::x(2),s(2)
  integer::i,n
  call rng_seed(rng,556677_i8)
  d=udmultivariate(2,bivnorm_pdf,bivnorm_logpdf)
  g=hitro_new(d,[0.0_dp,0.0_dp])
  do i=1,2000
    call g%sample(rng,x)
  end do
  s=0.0_dp
  n=12000
  do i=1,n
    call g%sample(rng,x)
    s=s+x
  end do
  s=s/real(n,dp)
  if(any(abs(s)>0.10_dp))then
    print *,'HITRO mean failed',s
    error stop 1
  end if
  if(g%acceptance_rate()<=0.05_dp .or. g%acceptance_rate()>=0.95_dp)then
    print *,'HITRO acceptance failed',g%acceptance_rate()
    error stop 1
  end if
  print *,'test_hitro: PASS'
end program test_hitro
