! SPDX-License-Identifier: GPL-3.0-only
module truncated_normal_linear
   use r_mod, only : dp
   use truncated_normal_math, only : lnNpr, var_tn
   implicit none
   private
   real(dp), parameter :: pi=acos(-1.0_dp)
   type, public :: cholperm_result
      real(dp), allocatable :: lmat(:,:), lower(:), upper(:)
      integer, allocatable :: perm(:)
      integer :: status=0
   end type
   public :: cholperm, dmvnorm, dmvt
contains

function lnNpr_scalar(a,b) result(w)
 real(dp),intent(in)::a,b
 real(dp)::w,t(1)
 t=lnNpr([a],[b])
 w=t(1)
end function

subroutine swap_int(a,b)
integer,intent(inout)::a,b
integer::t
t=a
a=b
b=t
end subroutine
subroutine swap_rows(a,i,j)
 real(dp),intent(inout)::a(:,:)
 integer,intent(in)::i,j
 real(dp)::t(size(a,2))
 if(i==j)return
 t=a(i,:)
 a(i,:)=a(j,:)
 a(j,:)=t
end subroutine

function trunc_mean(a,b,w) result(m)
 real(dp),intent(in)::a,b,w
 real(dp)::m,pa,pb
 pa=0.0_dp
 pb=0.0_dp
 if(abs(a)<huge(a)) pa=exp(-0.5_dp*a*a-w)/sqrt(2.0_dp*pi)
 if(abs(b)<huge(b)) pb=exp(-0.5_dp*b*b-w)/sqrt(2.0_dp*pi)
 m=pa-pb
end function

subroutine cholperm(sigma,l,u,res,method)
 real(dp),intent(in)::sigma(:,:),l(:),u(:)
 type(cholperm_result),intent(out)::res
 character(len=*),intent(in),optional::method
 integer::d,i,j,k,pos,idx,rel
 real(dp)::cii,mui,denom,w
 real(dp),allocatable::a(:),b(:),score(:),mu(:),ar(:),br(:),sr(:)
 character(len=3)::meth
 d=size(l)
 if(size(u)/=d.or.size(sigma,1)/=d.or.size(sigma,2)/=d) error stop 'cholperm: nonconformal input'
 meth='GGE'
 if(present(method))meth=method(1:min(3,len_trim(method)))
 allocate(res%lmat(d,d),res%lower(d),res%upper(d),res%perm(d),a(d),b(d),score(d),mu(d))
 res%lmat=0.0_dp
 mu=0.0_dp
 res%perm=[(i,i=1,d)]
 do i=1,d
  if(sigma(i,i)<=0.0_dp)then
  res%status=1
  return
  end if
  a(i)=l(i)/sqrt(sigma(i,i))
  b(i)=u(i)/sqrt(sigma(i,i))
 end do
 if(meth=='GB')then
 score=var_tn(a,b)
 else
 score=lnNpr(a,b)
 end if
 idx=minloc(score,dim=1)
 call swap_int(res%perm(1),res%perm(idx))
 cii=sqrt(sigma(res%perm(1),res%perm(1)))
 do i=1,d
 res%lmat(i,1)=sigma(res%perm(i),res%perm(1))/cii
 end do
 res%lmat(1,1)=cii
 w=lnNpr_scalar(a(res%perm(1)),b(res%perm(1)))
 mu(1)=trunc_mean(a(res%perm(1)),b(res%perm(1)),w)
 do j=2,d
  allocate(ar(d-j+1),br(d-j+1),sr(d-j+1))
  do i=j,d
   k=res%perm(i)
   mui=dot_product(res%lmat(i,1:j-1),mu(1:j-1))
   denom=sigma(k,k)-sum(res%lmat(i,1:j-1)**2)
   if(denom< -1e-8_dp)then
   res%status=2
   return
   end if
   denom=sqrt(max(denom,1e-20_dp))
   ar(i-j+1)=(l(k)-mui)/denom
   br(i-j+1)=(u(k)-mui)/denom
  end do
  if(meth=='GB')then
  sr=var_tn(ar,br)
  else
  sr=lnNpr(ar,br)
  end if
  rel=minloc(sr,dim=1)
  pos=j+rel-1
  if(pos/=j)then
  call swap_rows(res%lmat,j,pos)
  call swap_int(res%perm(j),res%perm(pos))
  end if
  k=res%perm(j)
  denom=sigma(k,k)-sum(res%lmat(j,1:j-1)**2)
  if(denom<=0.0_dp)then
  res%status=2
  return
  end if
  res%lmat(j,j)=sqrt(denom)
  do i=j+1,d
   idx=res%perm(i)
   res%lmat(i,j)=(sigma(idx,k)-dot_product(res%lmat(j,1:j-1),res%lmat(i,1:j-1)))/res%lmat(j,j)
  end do
  w=lnNpr_scalar(ar(rel),br(rel))
  mu(j)=trunc_mean(ar(rel),br(rel),w)
  deallocate(ar,br,sr)
 end do
 res%lower=l(res%perm)
 res%upper=u(res%perm)
end subroutine

subroutine cholesky_lower(a,l,info)
 real(dp),intent(in)::a(:,:)
 real(dp),intent(out)::l(:,:)
 integer,intent(out)::info
 integer::n,i,j
 n=size(a,1)
 l=0
 info=0
 do i=1,n
  do j=1,i
   if(i==j)then
   l(i,j)=a(i,i)-sum(l(i,1:j-1)**2)
   if(l(i,j)<=0)then
   info=i
   return
   end if
   l(i,j)=sqrt(l(i,j))
   else
   l(i,j)=(a(i,j)-dot_product(l(i,1:j-1),l(j,1:j-1)))/l(j,j)
   end if
  end do
 end do
end subroutine

subroutine forward_solve(l,b,x)
 real(dp),intent(in)::l(:,:),b(:)
 real(dp),intent(out)::x(:)
 integer::i
 do i=1,size(b)
 x(i)=(b(i)-dot_product(l(i,1:i-1),x(1:i-1)))/l(i,i)
 end do
end subroutine

function dmvnorm(x,mu,sigma,logd) result(out)
 real(dp),intent(in)::x(:,:),mu(:),sigma(:,:)
 logical,intent(in),optional::logd
 real(dp)::out(size(x,1)),l(size(mu),size(mu)),z(size(mu)),c
 integer::i,info
 logical::lg
 lg=.false.
 if(present(logd))lg=logd
 call cholesky_lower(sigma,l,info)
 if(info/=0)error stop 'dmvnorm: sigma not positive definite'
 c=-0.5_dp*size(mu)*log(2.0_dp*pi)-sum(log([(l(i,i),i=1,size(mu))]))
 do i=1,size(x,1)
 call forward_solve(l,x(i,:)-mu,z)
 out(i)=c-0.5_dp*sum(z*z)
 end do
 if(.not.lg)out=exp(out)
end function

function dmvt(x,mu,sigma,df,logd) result(out)
 real(dp),intent(in)::x(:,:),mu(:),sigma(:,:),df
 logical,intent(in),optional::logd
 real(dp)::out(size(x,1)),l(size(mu),size(mu)),z(size(mu)),c,p
 integer::i,info
 logical::lg
 lg=.false.
 if(present(logd))lg=logd
 p=real(size(mu),dp)
 call cholesky_lower(sigma,l,info)
 if(info/=0)error stop 'dmvt: sigma not positive definite'
 c=log_gamma(0.5_dp*(df+p))-log_gamma(0.5_dp*df)-0.5_dp*p*log(df*pi)-sum(log([(l(i,i),i=1,size(mu))]))
 do i=1,size(x,1)
 call forward_solve(l,x(i,:)-mu,z)
 out(i)=c-0.5_dp*(df+p)*log(1.0_dp+sum(z*z)/df)
 end do
 if(.not.lg)out=exp(out)
end function
end module
