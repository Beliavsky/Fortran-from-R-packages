program test_functions
  use trawl, only : dp,trawl_exp,trawl_dexp,trawl_supig,trawl_lm, &
    acf_exp,acf_dexp,acf_supig,acf_lm,lsd_mean,lsd_var,modlsd_mean, &
    modlsd_var,bivlsd_cor,bivlsd_cov,bivmodlsd_cov,bivmodlsd_cor
  implicit none
  integer :: fail
  real(dp) :: p1,p2,d,tp
  fail=0
  call chk(trawl_exp(-2.0_dp,0.5_dp),exp(-1.0_dp),1.0e-13_dp,'trawl_exp')
  call chk(trawl_dexp(-1.0_dp,0.3_dp,0.5_dp,2.0_dp), &
    0.3_dp*exp(-0.5_dp)+0.7_dp*exp(-2.0_dp),1.0e-13_dp,'trawl_dexp')
  call chk(acf_exp(10.0_dp,0.1_dp),exp(-1.0_dp),1.0e-13_dp,'acf_exp')
  call chk(acf_dexp(0.0_dp,0.1_dp,0.1_dp,1.0_dp),1.0_dp,1.0e-13_dp,'acf_dexp0')
  call chk(acf_supig(0.0_dp,0.1_dp,0.5_dp),1.0_dp,1.0e-13_dp,'acf_supig0')
  call chk(acf_lm(0.0_dp,0.1_dp,1.1_dp),1.0_dp,1.0e-13_dp,'acf_lm0')
  call chk(trawl_lm(-1.0_dp,0.5_dp,2.0_dp),1.0_dp/9.0_dp,1.0e-13_dp,'trawl_lm')
  call chk(trawl_supig(0.0_dp,0.2_dp,0.7_dp),1.0_dp,1.0e-13_dp,'trawl_supig0')
  call chk(lsd_mean(0.3_dp),0.3_dp/(0.7_dp*(-log(0.7_dp))),1.0e-13_dp,'lsd_mean')
  call assert_true(lsd_var(0.3_dp)>0.0_dp,'lsd_var positive')
  d=0.2_dp;tp=0.3_dp
  call chk(modlsd_mean(d,tp),(1.0_dp-d)*lsd_mean(tp),1.0e-13_dp,'mod mean')
  call assert_true(modlsd_var(d,tp)>0.0_dp,'mod var positive')
  p1=0.15_dp;p2=0.30_dp
  call assert_true(abs(bivlsd_cor(p1,p2))<1.0_dp,'biv cor bounds')
  call assert_true(abs(bivlsd_cov(p1,p2))>0.0_dp,'biv cov nonzero')
  call assert_true(abs(bivmodlsd_cov(0.2_dp,p1,p2))>0.0_dp,'biv mod cov nonzero')
  call assert_true(abs(bivmodlsd_cor(0.2_dp,p1,p2))<1.0_dp,'biv mod cor bounds')
  if(fail==0) then; print '(a)','test_functions: PASS'; else; error stop 1; end if
contains
  subroutine chk(a,b,tol,name)
    real(dp),intent(in)::a,b,tol;character(len=*),intent(in)::name
    if(abs(a-b)>tol) then; print *, 'FAIL ',trim(name),a,b;fail=fail+1;end if
  end subroutine
  subroutine assert_true(ok,name)
    logical,intent(in)::ok;character(len=*),intent(in)::name
    if(.not.ok) then;print *, 'FAIL ',trim(name);fail=fail+1;end if
  end subroutine
end program
