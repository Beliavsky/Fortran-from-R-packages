! SPDX-License-Identifier: GPL-3.0-only
program test_ordinary_approximations
  use poisson_binomial, only : dp, dpb_mean, dpb_geomean, dpb_geomean_counter, &
                               dpb_poisson, ppb_normal
  implicit none
  real(dp),parameter::p(4)=[0.1_dp,0.3_dp,0.6_dp,0.8_dp]
  real(dp),parameter::rmean(0:4)=[0.09150625_dp,0.299475_dp,0.3675375_dp, &
                                  0.200475_dp,0.04100625_dp]
  real(dp),parameter::rgeo(0:4)=[0.18248247641828605_dp,0.3868712786349386_dp, &
    0.3075693674201634_dp,0.10867687752661226_dp,0.0144_dp]
  real(dp),parameter::rgc(0:4)=[0.0504_dp,0.22388366713855232_dp, &
    0.37294565782296196_dp,0.2761125711534539_dp,0.07665810388503186_dp]
  real(dp),parameter::rpois(0:4)=[0.16529888822158653_dp,0.2975379987988558_dp, &
    0.26778419891897026_dp,0.16067051935138213_dp,0.10870839470920537_dp]
  real(dp),allocatable::a(:),cdf(:)
  a=dpb_mean(p); call chk(maxval(abs(a-rmean))<2e-13_dp,"mean")
  a=dpb_geomean(p); call chk(maxval(abs(a-rgeo))<2e-13_dp,"geomean")
  a=dpb_geomean_counter(p); call chk(maxval(abs(a-rgc))<2e-13_dp,"counter")
  a=dpb_poisson(p); call chk(maxval(abs(a-rpois))<3e-13_dp,"poisson")
  cdf=ppb_normal(p,.false.,.true.)
  call chk(abs(cdf(3)-0.7986081528767621_dp)<3e-13_dp,"normal cdf")
  cdf=ppb_normal(p,.true.,.true.)
  call chk(abs(cdf(3)-0.7988961652776708_dp)<3e-13_dp,"refined normal cdf")
  print '(a)', 'test_ordinary_approximations: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_ordinary_approximations
