! SPDX-License-Identifier: GPL-3.0-only
program test_var_t
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
  use imputefin
  use imputefin_rng, only : rng_state, rng_seed, rng_student_t
  implicit none
  integer,parameter::n=600,q=2
  real(dp)::y(n,q),ym(n,q),nanv,e1,e2
  type(rng_state)::rng
  type(var_t_result)::fit,fit_omit
  type(var_t_options)::opt
  integer::i
  nanv=ieee_value(0.0_dp,ieee_quiet_nan);call rng_seed(rng,8888_8)
  y(1,:)=[0.0_dp,0.0_dp]
  do i=2,n
    e1=0.18_dp*rng_student_t(rng,7.0_dp)
    e2=0.16_dp*(0.35_dp*e1/0.18_dp+sqrt(1.0_dp-0.35_dp**2)*rng_student_t(rng,7.0_dp))
    y(i,1)=0.03_dp+0.55_dp*y(i-1,1)+0.12_dp*y(i-1,2)+e1
    y(i,2)=-0.02_dp+0.08_dp*y(i-1,1)+0.45_dp*y(i-1,2)+e2
  end do
  ym=y;ym(100:106,1)=nanv;ym(240,2)=nanv;ym(400:404,:)=nanv
  opt%p=1;opt%maxiter=100;opt%tol=2.0e-4_dp
  call fit_var_t(ym,fit,opt)
  call check(fit%status==impute_ok.or.fit%status==impute_not_converged,'VAR fit status')
  call check(size(fit%phi,3)==1.and.size(fit%phi0)==q,'VAR result dimensions')
  call check(abs(fit%phi(1,1,1)-0.55_dp)<0.18_dp,'VAR phi11')
  call check(abs(fit%phi(2,2,1)-0.45_dp)<0.18_dp,'VAR phi22')
  call check(.not.any(ieee_is_nan(fit%completed)),'VAR completed data')
  opt%omit_missing=.true.
  call fit_var_t(ym,fit_omit,opt)
  call check(fit_omit%status==impute_ok.or.fit_omit%status==impute_not_converged,'VAR omit path')
  call check(fit_omit%n_used<n-1,'VAR omit rows')
  print '(a)', 'test_var_t: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//trim(msg);error stop 1;end if
  end subroutine check
end program test_var_t
