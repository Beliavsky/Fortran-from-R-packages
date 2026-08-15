module ld_random
use ld_kinds, only: dp, pi
use ld_linalg, only: chol_lower
implicit none
private
public :: seed_rng, rand_uniform, rand_normal, rand_gamma, rand_chisq, rand_beta
public :: rand_dirichlet, rand_mvn, rand_student_t, rand_mv_t, rand_categorical
contains
subroutine seed_rng(seed)
   integer, intent(in) :: seed
   integer :: n,i
   integer, allocatable :: put(:)
   call random_seed(size=n); allocate(put(n))
   do i=1,n; put(i)=mod(seed+104729*i,2147483646); if(put(i)<=0) put(i)=i; end do
   call random_seed(put=put)
end subroutine seed_rng

function rand_uniform() result(u)
   real(dp) :: u
   call random_number(u)
   u=max(u,tiny(1.0_dp)); u=min(u,1.0_dp-epsilon(1.0_dp))
end function rand_uniform

function rand_normal() result(z)
   real(dp) :: z,u1,u2
   u1=rand_uniform(); u2=rand_uniform()
   z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
end function rand_normal

recursive function rand_gamma(shape, scale) result(x)
   real(dp), intent(in) :: shape,scale
   real(dp) :: x,d,c,z,u,v
   if(shape<=0.0_dp .or. scale<=0.0_dp) then; x=0.0_dp; return; end if
   if(shape<1.0_dp) then
      x=rand_gamma(shape+1.0_dp,scale)*rand_uniform()**(1.0_dp/shape); return
   end if
   d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
   do
      z=rand_normal(); v=1.0_dp+c*z
      if(v<=0.0_dp) cycle
      v=v*v*v; u=rand_uniform()
      if(u<1.0_dp-0.0331_dp*z**4) exit
      if(log(u)<0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
   end do
   x=scale*d*v
end function rand_gamma

function rand_chisq(df) result(x)
   real(dp), intent(in) :: df
   real(dp) :: x
   x=rand_gamma(0.5_dp*df,2.0_dp)
end function rand_chisq

function rand_beta(a,b) result(x)
   real(dp), intent(in) :: a,b
   real(dp) :: x,g1,g2
   g1=rand_gamma(a,1.0_dp); g2=rand_gamma(b,1.0_dp); x=g1/(g1+g2)
end function rand_beta

subroutine rand_dirichlet(alpha,x)
   real(dp), intent(in) :: alpha(:)
   real(dp), intent(out) :: x(:)
   integer :: i
   do i=1,size(alpha); x(i)=rand_gamma(alpha(i),1.0_dp); end do
   x=x/sum(x)
end subroutine rand_dirichlet

subroutine rand_mvn(mu,sigma,x,info)
   real(dp), intent(in) :: mu(:),sigma(:,:)
   real(dp), intent(out) :: x(:)
   integer, intent(out), optional :: info
   real(dp), allocatable :: l(:,:),z(:)
   integer :: n,i,ifail
   n=size(mu); allocate(l(n,n),z(n)); call chol_lower(sigma,l,ifail)
   if(ifail/=0) then; x=mu; if(present(info)) info=ifail; return; end if
   do i=1,n; z(i)=rand_normal(); end do
   x=mu+matmul(l,z); if(present(info)) info=0
end subroutine rand_mvn

function rand_student_t(df) result(x)
   real(dp), intent(in) :: df
   real(dp) :: x
   x=rand_normal()/sqrt(rand_chisq(df)/df)
end function rand_student_t

subroutine rand_mv_t(mu,s,df,x,info)
   real(dp), intent(in) :: mu(:),s(:,:),df
   real(dp), intent(out) :: x(:)
   integer, intent(out), optional :: info
   real(dp) :: z(size(mu))
   integer :: ifail
   call rand_mvn(0.0_dp*mu,s,z,ifail)
   if(ifail==0) x=mu+z/sqrt(rand_chisq(df)/df); if(present(info)) info=ifail
end subroutine rand_mv_t

function rand_categorical(p) result(k)
   real(dp), intent(in) :: p(:)
   integer :: k,i
   real(dp) :: u,c,s
   s=sum(max(p,0.0_dp)); if(s<=0.0_dp) then; k=1; return; end if
   u=rand_uniform(); c=0.0_dp
   do i=1,size(p); c=c+max(p(i),0.0_dp)/s; if(u<=c) then; k=i; return; end if; end do
   k=size(p)
end function rand_categorical
end module ld_random
