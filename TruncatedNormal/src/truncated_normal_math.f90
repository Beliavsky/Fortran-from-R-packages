! Computational translation of TruncatedNormal 2.3.
! SPDX-License-Identifier: GPL-3.0-only
module truncated_normal_math
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_quiet_nan
   use r_compat, only : dp, normal_cdf, qnorm, runif1, dnorm
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: log_sqrt_2pi = 0.5_dp*log(2.0_dp*pi)
   public :: lnNpr, var_tn, phinv, norminvp, trandn, qtnorm_vec, rtnorm_vec
   public :: log_norm_cdf, log_norm_sf
contains

pure elemental function log1p_dp(x) result(v)
 real(dp),intent(in)::x
 real(dp)::v
 if(abs(x)<1.0e-5_dp)then
 v=x-x*x/2+x*x*x/3-x**4/4+x**5/5
 else
 v=log(1.0_dp+x)
 end if
end function
pure elemental function expm1_dp(x) result(v)
 real(dp),intent(in)::x
 real(dp)::v
 if(abs(x)<1.0e-5_dp)then
 v=x+x*x/2+x*x*x/6+x**4/24+x**5/120
 else
 v=exp(x)-1.0_dp
 end if
end function

pure elemental function log_norm_sf(x) result(v)
   real(dp), intent(in) :: x
   real(dp) :: v, xx, term
   if (x /= x) then
      v = x
   else if (x == ieee_value(x, ieee_positive_inf)) then
      v = ieee_value(x, ieee_negative_inf)
   else if (x == ieee_value(x, ieee_negative_inf)) then
      v = 0.0_dp
   else if (x < 8.0_dp) then
      v = log(0.5_dp*erfc(x/sqrt(2.0_dp)))
   else
      xx = x*x
      ! Mills-ratio expansion: Q(x)=phi(x)/x*(1-1/x^2+3/x^4-15/x^6+...)
      term = 1.0_dp - 1.0_dp/xx + 3.0_dp/(xx*xx) - 15.0_dp/(xx**3) + 105.0_dp/(xx**4)
      if (term <= 0.0_dp) term = 1.0_dp
      v = -0.5_dp*xx - log_sqrt_2pi - log(x) + log(term)
   end if
end function log_norm_sf

pure elemental function log_norm_cdf(x) result(v)
   real(dp), intent(in) :: x
   real(dp) :: v, q
   if (x < 0.0_dp) then
      v = log_norm_sf(-x)
   else if (x == ieee_value(x, ieee_positive_inf)) then
      v = 0.0_dp
   else
      q = 0.5_dp*erfc(x/sqrt(2.0_dp))
      v = log1p_dp(-q)
   end if
end function log_norm_cdf

pure elemental function logdiffexp(a, b) result(v)
   real(dp), intent(in) :: a, b
   real(dp) :: v
   ! log(exp(a)-exp(b))
   ! precondition a>=b.
   if (b == ieee_value(b, ieee_negative_inf)) then
      v = a
   else if (a == b) then
      v = ieee_value(a, ieee_negative_inf)
   else
      v = a + log1p_dp(-exp(b-a))
   end if
end function logdiffexp

pure function lnNpr(a, b) result(p)
   real(dp), intent(in) :: a(:), b(:)
   real(dp) :: p(size(a))
   integer :: i
   if (size(a) /= size(b)) error stop 'lnNpr: a and b must have equal length'
   if (any(a >= b)) error stop 'lnNpr: every lower bound must be less than upper bound'
   do i = 1, size(a)
      if (a(i) > 0.0_dp) then
         p(i) = logdiffexp(log_norm_sf(a(i)), log_norm_sf(b(i)))
      else if (b(i) < 0.0_dp) then
         p(i) = logdiffexp(log_norm_cdf(b(i)), log_norm_cdf(a(i)))
      else
         ! Crosses zero: 1 - Phi(a) - Q(b) = Phi(b)-Phi(a).
         p(i) = log1p_dp(-normal_cdf(a(i)) - 0.5_dp*erfc(b(i)/sqrt(2.0_dp)))
      end if
   end do
end function lnNpr

pure function var_tn(a, b) result(v)
   real(dp), intent(in) :: a(:), b(:)
   real(dp) :: v(size(a)), w(size(a)), pa, pb, mean_t, ta, tb
   integer :: i
   if (size(a) /= size(b)) error stop 'var_tn: a and b must have equal length'
   w = lnNpr(a,b)
   do i=1,size(a)
      if (ieee_is_finite(a(i))) then
         pa = exp(-0.5_dp*a(i)*a(i)-log_sqrt_2pi-w(i))
         ta = a(i)*pa
      else
         pa = 0.0_dp
         ta = 0.0_dp
      end if
      if (ieee_is_finite(b(i))) then
         pb = exp(-0.5_dp*b(i)*b(i)-log_sqrt_2pi-w(i))
         tb = b(i)*pb
      else
         pb = 0.0_dp
         tb = 0.0_dp
      end if
      mean_t = pa-pb
      v(i) = 1.0_dp + ta - tb - mean_t*mean_t
      if (v(i) < 0.0_dp .and. v(i) > -1.0e-12_dp) v(i)=0.0_dp
   end do
end function var_tn

pure elemental function phinv_scalar(p, l, u) result(x)
   real(dp), intent(in) :: p,l,u
   real(dp) :: x, ll, uu, pl, pu
   logical :: flip
   if (p <= 0.0_dp) then
      x=l
      return
   else if (p >= 1.0_dp) then
      x=u
      return
   end if
   flip = u < 0.0_dp
   if (flip) then
      ll = -u
      uu = -l
   else
      ll = l
      uu = u
   end if
   ! This path is used only away from the deep same-sign tails.
   pl = normal_cdf(ll)
   pu = normal_cdf(uu)
   x = qnorm(pl + (pu-pl)*p)
   if (flip) x=-x
end function phinv_scalar

pure function phinv(p,l,u) result(x)
   real(dp), intent(in) :: p(:),l(:),u(:)
   real(dp) :: x(size(p))
   integer :: i
   if (size(p)/=size(l) .or. size(p)/=size(u)) error stop 'phinv: nonconformal inputs'
   do i=1,size(p)
   x(i)=phinv_scalar(p(i),l(i),u(i))
   end do
end function phinv

pure elemental function qfun(x) result(q)
   real(dp), intent(in) :: x
   real(dp) :: q
   q = exp(0.5_dp*x*x + log_norm_sf(x))
end function qfun

pure elemental function normq_scalar(p,l,u) result(x)
   real(dp), intent(in) :: p,l,u
   real(dp) :: x, ql, qu, del, err, lsq, usq, val
   integer :: it
   if (l > 1.0e5_dp) then
      if (ieee_is_finite(u)) then
         x = sqrt(l*l - 2.0_dp*log1p_dp(p*expm1_dp(0.5_dp*(l*l-u*u))))
      else
         x = sqrt(l*l - 2.0_dp*log(1.0_dp-p))
      end if
      return
   end if
   ql=qfun(l)
   qu=0.0_dp
   if (ieee_is_finite(u)) qu=qfun(u)
   lsq=l*l
   if (ieee_is_finite(u)) then
      usq=u*u
      val=1.0_dp+p*expm1_dp(0.5_dp*(lsq-usq))
   else
      usq=huge(1.0_dp)
      val=1.0_dp-p
   end if
   val=max(val,tiny(1.0_dp))
   x=sqrt(max(0.0_dp,lsq-2.0_dp*log(val)))
   do it=1,100
      if (ieee_is_finite(u)) then
         del=-qfun(x)+(1.0_dp-p)*exp(0.5_dp*(x*x-lsq))*ql + p*exp(0.5_dp*(x*x-usq))*qu
      else
         del=-qfun(x)+(1.0_dp-p)*exp(0.5_dp*(x*x-lsq))*ql
      end if
      x=x-del
      err=abs(del)
      if (err<=1.0e-10_dp) exit
   end do
end function normq_scalar

pure elemental function norminvp_scalar(p,l,u) result(x)
   real(dp), intent(in) :: p,l,u
   real(dp) :: x
   if (p<=0.0_dp) then
      x=l
   else if (p>=1.0_dp) then
      x=u
   else if (l>35.0_dp) then
      x=normq_scalar(p,l,u)
   else if (u< -35.0_dp) then
      x=-normq_scalar(1.0_dp-p,-u,-l)
   else
      x=phinv_scalar(p,l,u)
   end if
end function norminvp_scalar

pure function norminvp(p,l,u) result(x)
   real(dp), intent(in) :: p(:),l(:),u(:)
   real(dp) :: x(size(p))
   integer :: i
   if(size(p)/=size(l).or.size(p)/=size(u)) error stop 'norminvp: nonconformal inputs'
   if(any(l>u).or.any(p<0.0_dp).or.any(p>1.0_dp)) error stop 'norminvp: invalid bounds/probabilities'
   do i=1,size(p)
   x(i)=norminvp_scalar(p(i),l(i),u(i))
   end do
end function norminvp

function ntail_scalar(l,u) result(x)
   real(dp), intent(in) :: l,u
   real(dp) :: x,c,f,y
   c=0.5_dp*l*l
   if (ieee_is_finite(u)) then
      f=expm1_dp(c-0.5_dp*u*u)
   else
      f=-1.0_dp
   end if
   do
      y=c-log1p_dp(runif1()*f)
      if(runif1()*runif1()*y <= c) exit
   end do
   x=sqrt(2.0_dp*y)
end function ntail_scalar

function trnd_scalar(l,u) result(x)
   real(dp), intent(in) :: l,u
   real(dp) :: x
   do
      x = qnorm(runif1())
      if (x>=l .and. x<=u) exit
   end do
end function trnd_scalar

function tn_scalar(l,u) result(x)
   real(dp), intent(in) :: l,u
   real(dp) :: x,pl,pu
   if(abs(u-l)>2.05_dp) then
      x=trnd_scalar(l,u)
   else
      pl=normal_cdf(l)
      pu=normal_cdf(u)
      x=qnorm(pl+(pu-pl)*runif1())
   end if
end function tn_scalar

function trandn(l,u) result(x)
   real(dp), intent(in) :: l(:),u(:)
   real(dp) :: x(size(l))
   integer :: i
   if(size(l)/=size(u).or.any(l>u)) error stop 'trandn: invalid bounds'
   do i=1,size(l)
      if(l(i)>0.4_dp) then
         x(i)=ntail_scalar(l(i),u(i))
      else if(u(i)<-0.4_dp) then
         x(i)=-ntail_scalar(-u(i),-l(i))
      else
         x(i)=tn_scalar(l(i),u(i))
      end if
   end do
end function trandn

pure function qtnorm_vec(p,mu,sd,lb,ub) result(x)
   real(dp),intent(in)::p(:),mu(:),sd(:),lb(:),ub(:)
   real(dp)::x(size(p)),l2(size(p)),u2(size(p))
   if(any([size(mu)/=size(p),size(sd)/=size(p),size(lb)/=size(p),size(ub)/=size(p)])) error stop 'qtnorm_vec: sizes differ'
   if(any(sd<=0.0_dp).or.any(lb>=ub)) error stop 'qtnorm_vec: invalid sd/bounds'
   l2=(lb-mu)/sd
   u2=(ub-mu)/sd
   x=mu+sd*norminvp(p,l2,u2)
end function qtnorm_vec

function rtnorm_vec(n,mu,sd,lb,ub,fast) result(x)
   integer,intent(in)::n
   real(dp),intent(in)::mu(:),sd(:),lb(:),ub(:)
   logical,intent(in),optional::fast
   real(dp),allocatable::x(:,:)
   real(dp),allocatable::l2(:),u2(:),z(:),p(:)
   integer::d,i,j
   logical::use_fast
   d=size(mu)
   if(size(sd)/=d.or.size(lb)/=d.or.size(ub)/=d) error stop 'rtnorm_vec: sizes differ'
   if(any(sd<=0.0_dp).or.any(lb>=ub).or.n<1) error stop 'rtnorm_vec: invalid input'
   use_fast=.true.
   if(present(fast)) use_fast=fast
   allocate(x(n,d),l2(d),u2(d),z(d),p(d))
   l2=(lb-mu)/sd
   u2=(ub-mu)/sd
   do i=1,n
      if(use_fast) then
         z=trandn(l2,u2)
      else
         do j=1,d
         p(j)=runif1()
         end do
         z=norminvp(p,l2,u2)
      end if
      x(i,:)=mu+sd*z
   end do
end function rtnorm_vec

end module truncated_normal_math
