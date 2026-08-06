! SPDX-License-Identifier: GPL-3.0-only
program test_classical_alpha
  use spantest
  implicit none
  type(simulation_result) :: sim
  type(span_result) :: bj,grs,f1,py,bj_alt,grs_alt
  real(dp), allocatable :: alt(:,:),bad(:,:)

  sim=span_simulate(250,3,5,ncp=0.0_dp,dgp=1,seed=42)
  call check(sim%status==span_ok,'simulation failed')
  bj=span_bj(sim%r1,sim%r2)
  grs=span_grs(sim%r1,sim%r2)
  f1=span_f1(sim%r1,sim%r2)
  py=span_py(sim%r1,sim%r2)
  call check(bj%status==span_ok .and. grs%status==span_ok,'BJ/GRS failed')
  call check(f1%status==span_ok .and. py%status==span_ok,'F1/PY failed')
  call check(abs(bj%stat-grs%stat)<1.0e-10_dp,'BJ and GRS statistics differ')
  call check(abs(bj%pval-grs%pval)<1.0e-10_dp,'BJ and GRS p-values differ')
  call check(in_unit(bj%pval) .and. in_unit(f1%pval) .and. in_unit(py%pval),'invalid p-value')

  allocate(alt(size(sim%r2,1),size(sim%r2,2)))
  alt=sim%r2+0.8_dp
  bj_alt=span_bj(sim%r1,alt)
  grs_alt=span_grs(sim%r1,alt)
  call check(bj_alt%pval<1.0e-8_dp .and. grs_alt%pval<1.0e-8_dp,'alpha alternative not detected')

  allocate(bad(size(sim%r1,1),2))
  bad(:,1)=sim%r1(:,1); bad(:,2)=sim%r1(:,1)
  bj=span_bj(bad,sim%r2)
  call check(bj%status==span_singular,'singular design not detected')
  print '(a)','test_classical_alpha: PASS'
contains
  pure logical function in_unit(x)
    real(dp),intent(in)::x
    in_unit=x>=0.0_dp .and. x<=1.0_dp
  end function
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if (.not.ok) error stop msg
  end subroutine
end program test_classical_alpha
