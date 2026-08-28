! SPDX-License-Identifier: GPL-3.0-only
module truncated_normal_tilting
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use r_mod, only : dp, dnorm
   use truncated_normal_math, only : lnNpr
   implicit none
   private
   real(dp),parameter::pi=acos(-1.0_dp)
   public :: gradpsi_eval,jacpsi_eval,psy_eval,gradpsi_t_eval,psy_t_eval
contains
function lnNpr_scalar(a,b) result(w)
 real(dp),intent(in)::a,b
 real(dp)::w,t(1)
 t=lnNpr([a],[b])
 w=t(1)
end function
subroutine gradpsi_eval(y,lmat,l,u,g)
 real(dp),intent(in)::y(:),lmat(:,:),l(:),u(:)
 real(dp),intent(out)::g(:)
 integer::d
 real(dp)::x(size(l)),mu(size(l)),cv(size(l)),lt(size(l)),ut(size(l)),w(size(l)),pl(size(l)),pu(size(l)),P(size(l))
 d=size(l)
 x=0
 mu=0
 x(1:d-1)=y(1:d-1)
 mu(1:d-1)=y(d:2*d-2)
 cv=matmul(lmat,x)
 lt=l-mu-cv
 ut=u-mu-cv
 w=lnNpr(lt,ut)
 pl=0
 pu=0
 where(ieee_is_finite(lt))pl=exp(-0.5_dp*lt*lt-w)/sqrt(2.0_dp*pi)
 where(ieee_is_finite(ut))pu=exp(-0.5_dp*ut*ut-w)/sqrt(2.0_dp*pi)
 P=pl-pu
 g(1:d-1)=-mu(1:d-1)+matmul(transpose(lmat(:,1:d-1)),P)
 g(d:2*d-2)=mu(1:d-1)-x(1:d-1)+P(1:d-1)
end subroutine
subroutine jacpsi_eval(y,lmat,l,u,J)
 real(dp),intent(in)::y(:),lmat(:,:),l(:),u(:)
 real(dp),intent(out)::J(:,:)
 integer::d,i
 real(dp) :: x(size(l)), mu(size(l)), cv(size(l)), lt(size(l)), ut(size(l))
 real(dp) :: w(size(l)), pl(size(l)), pu(size(l)), P(size(l)), dpv(size(l))
 real(dp) :: DL(size(lmat,1),size(lmat,2)), mx(size(lmat,1),size(lmat,2))
 real(dp) :: xx(size(lmat,2),size(lmat,2))
 d=size(l)
 x=0
 mu=0
 x(1:d-1)=y(1:d-1)
 mu(1:d-1)=y(d:2*d-2)
 cv=matmul(lmat,x)
 lt=l-mu-cv
 ut=u-mu-cv
 w=lnNpr(lt,ut)
 pl=0
 pu=0
 where(ieee_is_finite(lt))pl=exp(-0.5_dp*lt*lt-w)/sqrt(2.0_dp*pi)
 where(ieee_is_finite(ut))pu=exp(-0.5_dp*ut*ut-w)/sqrt(2.0_dp*pi)
 P=pl-pu
 where(.not.ieee_is_finite(lt))lt=0
 where(.not.ieee_is_finite(ut))ut=0
 dpv=-P*P+lt*pl-ut*pu
 do i=1,d
 DL(i,:)=dpv(i)*lmat(i,:)
 end do
 mx=DL
 do i=1,d
 mx(i,i)=mx(i,i)-1
 end do
 xx=matmul(transpose(lmat),DL)
 J(1:d-1,1:d-1)=xx(1:d-1,1:d-1)
 J(1:d-1,d:2*d-2)=transpose(mx(1:d-1,1:d-1))
 J(d:2*d-2,1:d-1)=mx(1:d-1,1:d-1)
 J(d:2*d-2,d:2*d-2)=0
 do i=1,d-1
 J(d+i-1,d+i-1)=1+dpv(i)
 end do
end subroutine
function psy_eval(xin,lmat,l,u,muin) result(v)
 real(dp),intent(in)::xin(:),lmat(:,:),l(:),u(:),muin(:)
 real(dp)::v,x(size(l)),mu(size(l)),cv(size(l))
 x=0
 mu=0
 x(1:min(size(xin),size(l)-1))=xin(1:min(size(xin),size(l)-1))
 mu(1:min(size(muin),size(l)-1))=muin(1:min(size(muin),size(l)-1))
 cv=matmul(lmat,x)
 v=sum(lnNpr(l-mu-cv,u-mu-cv)+0.5_dp*mu*mu-x*mu)
end function
subroutine gradpsi_t_eval(y,lmat,l,u,nu,g)
 real(dp),intent(in)::y(:),lmat(:,:),l(:),u(:),nu
 real(dp),intent(out)::g(:)
 integer::d
 real(dp) :: x(size(l)), mu(size(l)), co(size(l)), lt(size(l)), ut(size(l))
 real(dp) :: w(size(l)), pl(size(l)), pu(size(l)), P(size(l)), ls(size(l)), us(size(l))
 real(dp) :: r, eta
 d=size(l)
 x=0
 mu=0
 x(1:d-1)=y(1:d-1)
 r=exp(y(d))
 mu(1:d-1)=y(d+1:2*d-1)
 eta=y(2*d)
 ls=l/sqrt(nu)
 us=u/sqrt(nu)
 co=matmul(lmat,x)
 lt=r*ls-mu-co
 ut=r*us-mu-co
 w=lnNpr(lt,ut)
 pl=0
 pu=0
 where(ieee_is_finite(lt))pl=exp(-0.5_dp*lt*lt-w)/sqrt(2.0_dp*pi)
 where(ieee_is_finite(ut))pu=exp(-0.5_dp*ut*ut-w)/sqrt(2.0_dp*pi)
 P=pl-pu
 g(1:d-1)=-mu(1:d-1)+matmul(transpose(lmat(:,1:d-1)),P)
 g(d)= (nu-1)/r-eta+sum(merge(us*pu,0.0_dp,ieee_is_finite(us))-merge(ls*pl,0.0_dp,ieee_is_finite(ls)))
 g(d+1:2*d-1)=mu(1:d-1)-x(1:d-1)+P(1:d-1)
 g(2*d)=eta-r+exp(dnorm(eta,log_=.true.)-lnNpr_scalar(-eta,huge(eta)))
end subroutine
function psy_t_eval(xin,lmat,l,u,nu,muin) result(v)
 real(dp),intent(in)::xin(:),lmat(:,:),l(:),u(:),nu,muin(:)
 real(dp)::v,x(size(l)),mu(size(l)),cv(size(l)),r,eta,w1(1)
 integer::d
 d=size(l)
 x=0
 mu=0
 x(1:d)=xin(1:d)
 mu(1:d)=muin(1:d)
 r=x(d)
 eta=mu(d)
 x(d)=0
 mu(d)=0
 cv=matmul(lmat,x)
 w1=lnNpr([-eta],[huge(eta)])
 v = sum(lnNpr(r*l/sqrt(nu)-mu-cv, r*u/sqrt(nu)-mu-cv)) &
   + 0.5_dp*sum(mu*mu) - sum(x*mu) + 0.5_dp*log(2*pi) - log_gamma(nu/2) &
   - (nu/2-1)*log(2.0_dp) + 0.5_dp*eta*eta - r*eta + (nu-1)*log(r) + w1(1)
end function
end module
