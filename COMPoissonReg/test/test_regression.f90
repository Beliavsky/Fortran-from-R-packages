program test_regression
   use compoissonreg
   implicit none
   integer,parameter :: n=10,nz=20
   integer :: y(n), yz(nz), fails
   real(dp) :: x(n,1),s(n,1),xz(nz,1),sz(nz,1),wz(nz,1)
   type(cmp_fit_t) :: fit
   type(zicmp_fit_t) :: zfit
   type(cmp_init_t) :: ini
   type(cmp_fixed_t) :: fix
   real(dp) :: target_beta, p_hat, target_zeta
   fails=0
   y=[0,1,2,1,3,2,4,1,0,2]
   x=1.0_dp;s=1.0_dp
   ini=default_init(1,1);fix=default_fixed(1,1);fix%gamma=.true.
   call fit_cmp_raw(y,x,s,fit,ini,fix)
   target_beta=log(sum(real(y,dp))/real(n,dp))
   call check_close('Poisson regression beta',fit%beta(1),target_beta,2.0e-4_dp,fails)
   if(.not.fit%converged)then;print *,'FAIL CMP convergence';fails=fails+1;end if

   yz=[0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,3,3,4,1,2]
   xz=1.0_dp;sz=1.0_dp;wz=1.0_dp
   ini=default_init(1,1,1);ini%beta(1)=log(2.0_dp);ini%gamma(1)=0.0_dp;ini%zeta(1)=0.0_dp
   fix=default_fixed(1,1,1);fix%beta=.true.;fix%gamma=.true.
   call fit_zicmp_raw(yz,xz,sz,wz,zfit,ini,fix)
   p_hat=(0.4_dp-exp(-2.0_dp))/(1.0_dp-exp(-2.0_dp))
   target_zeta=log(p_hat/(1.0_dp-p_hat))
   call check_close('ZICMP zeta',zfit%zeta(1),target_zeta,5.0e-3_dp,fails)
   if(.not.zfit%converged)then;print *,'FAIL ZICMP convergence';fails=fails+1;end if

   if(fails==0)then;print *,'test_regression: PASS';else;error stop 1;end if
contains
   subroutine check_close(name,a,b,tol,fails)
      character(*),intent(in)::name;real(dp),intent(in)::a,b,tol;integer,intent(inout)::fails
      if(abs(a-b)>tol*max(1.0_dp,abs(b)))then;print *,'FAIL ',trim(name),a,b;fails=fails+1;end if
   end subroutine check_close
end program test_regression
