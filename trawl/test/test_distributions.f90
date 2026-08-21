program test_distributions
  use trawl, only : dp,set_trawl_seed,bivariate_nbsim,bivariate_lsdsim,trivariate_lsdsim, &
    modlsd_mean,modlsd_var,bivlsd_cor
  implicit none
  integer,allocatable :: x(:,:)
  real(dp),allocatable :: a(:),b(:)
  real(dp)::p1,p2,mu1,mu2,v1,v2,corr
  integer::n,s,fail
  fail=0;call set_trawl_seed(13579)
  n=30000;p1=0.1_dp;p2=0.85_dp
  call bivariate_nbsim(n,3.0_dp,p1,p2,x,s)
  call assert_true(s==0 .and. size(x,1)==n,'bivariate NB status')
  allocate(a(n),b(n));a=real(x(:,1),dp);b=real(x(:,2),dp)
  mu1=3.0_dp*p1/(1.0_dp-p1-p2);mu2=3.0_dp*p2/(1.0_dp-p1-p2)
  v1=3.0_dp*p1*(1.0_dp-p2)/(1.0_dp-p1-p2)**2
  v2=3.0_dp*p2*(1.0_dp-p1)/(1.0_dp-p1-p2)**2
  call relchk(mean(a),mu1,0.04_dp,'NB mean1');call relchk(mean(b),mu2,0.025_dp,'NB mean2')
  call relchk(var(a),v1,0.08_dp,'NB var1');call relchk(var(b),v2,0.06_dp,'NB var2')
  deallocate(x,a,b)

  n=12000;p1=0.15_dp;p2=0.30_dp
  call bivariate_lsdsim(n,p1,p2,x,s);allocate(a(n),b(n));a=real(x(:,1),dp);b=real(x(:,2),dp)
  mu1=modlsd_mean(log(1.0_dp-p2)/log(1.0_dp-p1-p2),p1/(1.0_dp-p2))
  mu2=modlsd_mean(log(1.0_dp-p1)/log(1.0_dp-p1-p2),p2/(1.0_dp-p1))
  call abscheck(mean(a),mu1,0.06_dp,'LSD mean1');call abscheck(mean(b),mu2,0.06_dp,'LSD mean2')
  corr=cor(a,b);call abscheck(corr,bivlsd_cor(p1,p2),0.08_dp,'LSD cor')
  deallocate(x,a,b)

  call trivariate_lsdsim(3000,0.15_dp,0.25_dp,0.55_dp,x,s)
  call assert_true(s==0 .and. size(x,2)==3,'trivariate LSD status')
  call assert_true(minval(x)>=0,'trivariate nonnegative')
  if(fail==0) then;print '(a)','test_distributions: PASS';else;error stop 1;end if
contains
  real(dp) function mean(z) result(m);real(dp),intent(in)::z(:);m=sum(z)/size(z);end function
  real(dp) function var(z) result(v);real(dp),intent(in)::z(:);real(dp)::m;m=mean(z);v=sum((z-m)**2)/(size(z)-1);end function
  real(dp) function cor(z,w) result(r);real(dp),intent(in)::z(:),w(:);real(dp)::mz,mw
    mz=mean(z);mw=mean(w);r=sum((z-mz)*(w-mw))/sqrt(sum((z-mz)**2)*sum((w-mw)**2));end function
  subroutine relchk(a0,b0,tol,name);real(dp),intent(in)::a0,b0,tol;character(len=*),intent(in)::name
    if(abs(a0-b0)>tol*max(1.0_dp,abs(b0)))then;print *, 'FAIL ',trim(name),a0,b0;fail=fail+1;end if;end subroutine
  subroutine abscheck(a0,b0,tol,name);real(dp),intent(in)::a0,b0,tol;character(len=*),intent(in)::name
    if(abs(a0-b0)>tol)then;print *, 'FAIL ',trim(name),a0,b0;fail=fail+1;end if;end subroutine
  subroutine assert_true(ok,name);logical,intent(in)::ok;character(len=*),intent(in)::name
    if(.not.ok)then;print *, 'FAIL ',trim(name);fail=fail+1;end if;end subroutine
end program
