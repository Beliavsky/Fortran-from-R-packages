module aep_distribution
   use aep_special, only: dp, reg_gamma_p, gamma_quantile
   implicit none
   private
   public :: daep, paep, qaep, raep
   interface daep; module procedure daep_scalar, daep_vec; end interface
   interface paep; module procedure paep_scalar, paep_vec; end interface
   interface qaep; module procedure qaep_scalar, qaep_vec; end interface
contains
   pure real(dp) function daep_scalar(x, alpha, sigma, mu, epsilon, log_pdf) result(v)
      real(dp), intent(in) :: x,alpha,sigma,mu,epsilon
      logical, intent(in), optional :: log_pdf
      logical :: lp
      real(dp) :: den, z
      lp=.false.; if (present(log_pdf)) lp=log_pdf
      if (.not.valid_par(alpha,sigma,epsilon)) then; v=nanv(); return; end if
      den = sigma*(1.0_dp + sign(1.0_dp,x-mu)*epsilon)
      if (x==mu) den=sigma*(1.0_dp+epsilon)
      z=abs(x-mu)/den
      v=-log(2.0_dp*sigma)-log_gamma(1.0_dp+1.0_dp/alpha)-z**alpha
      if (.not.lp) v=exp(v)
   end function
   pure function daep_vec(x,alpha,sigma,mu,epsilon,log_pdf) result(v)
      real(dp), intent(in)::x(:),alpha,sigma,mu,epsilon
      logical,intent(in),optional::log_pdf
      real(dp)::v(size(x)); integer::i
      do i=1,size(x); v(i)=daep_scalar(x(i),alpha,sigma,mu,epsilon,log_pdf); end do
   end function

   pure real(dp) function paep_scalar(x,alpha,sigma,mu,epsilon,log_p,lower_tail) result(v)
      real(dp),intent(in)::x,alpha,sigma,mu,epsilon
      logical,intent(in),optional::log_p,lower_tail
      logical::lp,lt; real(dp)::g
      lp=.false.; lt=.true.; if(present(log_p))lp=log_p; if(present(lower_tail))lt=lower_tail
      if(.not.valid_par(alpha,sigma,epsilon))then;v=nanv();return;end if
      if(x<=mu)then
         g=reg_gamma_p(1.0_dp/alpha,((mu-x)/(sigma*(1.0_dp-epsilon)))**alpha)
         v=0.5_dp*(1.0_dp-epsilon)*(1.0_dp-g)
      else
         g=reg_gamma_p(1.0_dp/alpha,((x-mu)/(sigma*(1.0_dp+epsilon)))**alpha)
         v=0.5_dp*(1.0_dp-epsilon)+0.5_dp*(1.0_dp+epsilon)*g
      end if
      if(.not.lt)v=1.0_dp-v
      if(lp)v=log(v)
   end function
   pure function paep_vec(x,alpha,sigma,mu,epsilon,log_p,lower_tail) result(v)
      real(dp),intent(in)::x(:),alpha,sigma,mu,epsilon; logical,intent(in),optional::log_p,lower_tail
      real(dp)::v(size(x));integer::i
      do i=1,size(x);v(i)=paep_scalar(x(i),alpha,sigma,mu,epsilon,log_p,lower_tail);end do
   end function

   pure real(dp) function qaep_scalar(u,alpha,sigma,mu,epsilon) result(v)
      real(dp),intent(in)::u,alpha,sigma,mu,epsilon
      real(dp)::q,cut
      if(.not.valid_par(alpha,sigma,epsilon).or.u<0.0_dp.or.u>1.0_dp)then;v=nanv();return;end if
      if(u==0.0_dp)then;v=-huge(1.0_dp);return;else if(u==1.0_dp)then;v=huge(1.0_dp);return;end if
      cut=0.5_dp*(1.0_dp-epsilon)
      if(u<cut)then
         q=gamma_quantile(1.0_dp-2.0_dp*u/(1.0_dp-epsilon),1.0_dp/alpha)
         v=mu-sigma*(1.0_dp-epsilon)*q**(1.0_dp/alpha)
      else
         q=gamma_quantile(2.0_dp*(u-cut)/(1.0_dp+epsilon),1.0_dp/alpha)
         v=mu+sigma*(1.0_dp+epsilon)*q**(1.0_dp/alpha)
      end if
   end function
   pure function qaep_vec(u,alpha,sigma,mu,epsilon) result(v)
      real(dp),intent(in)::u(:),alpha,sigma,mu,epsilon;real(dp)::v(size(u));integer::i
      do i=1,size(u);v(i)=qaep_scalar(u(i),alpha,sigma,mu,epsilon);end do
   end function

   subroutine raep(x,alpha,sigma,mu,epsilon)
      real(dp),intent(out)::x(:);real(dp),intent(in)::alpha,sigma,mu,epsilon
      real(dp)::u(size(x));integer::i
      call random_number(u)
      do i=1,size(x);x(i)=qaep_scalar(u(i),alpha,sigma,mu,epsilon);end do
   end subroutine

   pure logical function valid_par(alpha,sigma,epsilon)
      real(dp),intent(in)::alpha,sigma,epsilon
      valid_par=alpha>0.0_dp.and.sigma>0.0_dp.and.abs(epsilon)<1.0_dp
   end function
   pure real(dp) function nanv() result(x)
      use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function
end module aep_distribution
