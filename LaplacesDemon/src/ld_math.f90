module ld_math
use ld_kinds, only: dp
implicit none
private
public :: logit, invlogit, cloglog, invcloglog, loglog, invloglog
public :: center_scale, precision_to_variance, variance_to_precision
contains
pure function logit(p) result(x)
   real(dp),intent(in)::p; real(dp)::x
   x=log(p/(1.0_dp-p))
end function logit
pure function invlogit(x) result(p)
   real(dp),intent(in)::x; real(dp)::p
   if(x>=0.0_dp) then; p=1.0_dp/(1.0_dp+exp(-x)); else; p=exp(x)/(1.0_dp+exp(x)); end if
end function invlogit
pure function cloglog(p) result(x)
   real(dp),intent(in)::p; real(dp)::x
   x=log(-log(1.0_dp-p))
end function cloglog
pure function invcloglog(x) result(p)
   real(dp),intent(in)::x; real(dp)::p
   p=1.0_dp-exp(-exp(x))
end function invcloglog
pure function loglog(p) result(x)
   real(dp),intent(in)::p; real(dp)::x
   x=-log(-log(p))
end function loglog
pure function invloglog(x) result(p)
   real(dp),intent(in)::x; real(dp)::p
   p=exp(-exp(-x))
end function invloglog
subroutine center_scale(x,z,means,sds)
   real(dp),intent(in)::x(:,:); real(dp),intent(out)::z(:,:),means(:),sds(:)
   integer::j,n
   real(dp)::m
   n=size(x,1)
   do j=1,size(x,2)
      m=sum(x(:,j))/real(n,dp); means(j)=m
      if(n>1) then; sds(j)=sqrt(sum((x(:,j)-m)**2)/real(n-1,dp)); else; sds(j)=1.0_dp; end if
      if(sds(j)<=0.0_dp) sds(j)=1.0_dp; z(:,j)=(x(:,j)-m)/sds(j)
   end do
end subroutine center_scale
pure function precision_to_variance(tau) result(v)
   real(dp),intent(in)::tau; real(dp)::v; v=1.0_dp/tau
end function precision_to_variance
pure function variance_to_precision(v) result(tau)
   real(dp),intent(in)::v; real(dp)::tau; tau=1.0_dp/v
end function variance_to_precision
end module ld_math
