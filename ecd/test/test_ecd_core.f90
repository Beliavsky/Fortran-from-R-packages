! SPDX-License-Identifier: Artistic-2.0
program test_ecd_core
  use ecd_api
  implicit none
  type(ecd_model) :: d
  type(ecld_model) :: l
  type(ecd_stats_type) :: s
  real(dp) :: q
  integer :: st

  d=ecd_new(with_stats=.true.,status=st)
  call check(st==ecd_ok,'standard cusp construction')
  call close(d%norm_const,1.5_dp*sqrt_pi,2.0e-12_dp,'standard cusp constant')
  call close(ecd_pdf(d,0.0_dp),1.0_dp/(1.5_dp*sqrt_pi),2.0e-12_dp,'standard cusp pdf')
  call close(ecd_cdf(d,0.0_dp),0.5_dp,5.0e-10_dp,'standard cusp cdf')
  call close(ecd_cusp_std_moment(2),13.125_dp,2.0e-12_dp,'standard cusp second moment')
  call close(ecd_cusp_std_moment(4)/ecd_cusp_std_moment(2)**2,429.0_dp/35.0_dp, &
    2.0e-12_dp,'standard cusp kurtosis')
  q=ecd_quantile(d,0.9_dp,st)
  call check(st==ecd_ok,'standard cusp quantile status')
  call close(ecd_cdf(d,q),0.9_dp,5.0e-9_dp,'standard cusp inversion')
  s=ecd_statistics(d)
  call close(s%variance,13.125_dp,5.0e-9_dp,'standard cusp numerical variance')
  call close(real(ecd_cusp_std_cf(0.0_dp),dp),1.0_dp,1.0e-15_dp,'standard cusp cf zero')

  call close(ecd_cusp_a2r(2.0_dp),-(27.0_dp)**(1.0_dp/3.0_dp),2.0e-14_dp,'cusp alpha to gamma')
  call close(ecd_cusp_r2a(ecd_cusp_a2r(2.0_dp)),2.0_dp,2.0e-13_dp,'cusp conversion inverse')

  l=ecld_new(lambda=3.0_dp,sigma=0.4_dp)
  call close(ecld_const(l),1.06347231054331_dp,2.0e-13_dp,'ecld constant')
  call close(ecld_pdf(l,0.2_dp),0.5008243470289347_dp,3.0e-12_dp,'ecld pdf reference')
  call close(ecld_cdf(l,0.2_dp),0.6306642537556907_dp,3.0e-9_dp,'ecld cdf reference')
  call close(ecld_variance(l),2.1_dp,2.0e-13_dp,'ecld variance')
  call close(ecld_kurtosis(l),12.257142857142857_dp,2.0e-12_dp,'ecld kurtosis')
  call close(ecld_mgf(ecld_new(lambda=1.0_dp,sigma=0.4_dp),0.7_dp), &
    exp(0.4_dp**2*0.7_dp**2/4.0_dp),2.0e-12_dp,'gaussian ecld mgf')

  print '(a)', 'test_ecd_core: PASS'
contains
  subroutine close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    if(abs(x-y)>tol*max(1.0_dp,abs(y)))then
      write(*,*)trim(msg),x,y; error stop 1
    end if
  end subroutine close
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,*)trim(msg);error stop 1;end if
  end subroutine check
end program test_ecd_core
