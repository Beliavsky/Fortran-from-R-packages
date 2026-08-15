module ld_quadrature_catalog
use ld_kinds, only: dp
use ld_interfaces, only: log_target_iface
use ld_linalg, only: chol_lower, make_positive_definite, outer_product
use ld_quadrature, only: gauss_hermite_rule
use ld_evidence, only: iterative_quad_result_t
implicit none
private
public :: componentwise_iterative_quadrature, adaptive_sparse_grid_quadrature

contains

subroutine componentwise_iterative_quadrature(f,mean0,var0,order,max_iter,tol,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: mean0(:),var0(:)
   integer, intent(in) :: order,max_iter
   real(dp), intent(in) :: tol
   type(iterative_quad_result_t), intent(out) :: res
   integer :: p,it,j,k
   real(dp), allocatable :: nodes(:),ghw(:),lw(:),w(:),theta(:),newm(:),newv(:)
   real(dp) :: mlog,lse,scale,change,den
   p=size(mean0); allocate(nodes(order),ghw(order),lw(order),w(order),theta(p),newm(p),newv(p))
   allocate(res%mean(p),res%covariance(p,p)); call gauss_hermite_rule(order,nodes,ghw)
   res%mean=mean0; res%covariance=0.0_dp
   do j=1,p; res%covariance(j,j)=max(var0(j),1.0e-10_dp); end do
   res%log_normalizer=0.0_dp
   do it=1,max_iter
      newm=res%mean; newv=[(res%covariance(j,j),j=1,p)]; res%log_normalizer=0.0_dp
      do j=1,p
         scale=sqrt(max(res%covariance(j,j),1.0e-12_dp)); theta=res%mean
         do k=1,order
            theta(j)=res%mean(j)+sqrt(2.0_dp)*scale*nodes(k)
            lw(k)=log(max(ghw(k),tiny(1.0_dp)))+f(theta)+nodes(k)*nodes(k)
         end do
         mlog=maxval(lw); w=exp(lw-mlog); den=sum(w); if(den<=0.0_dp) cycle; w=w/den
         newm(j)=sum(w*(res%mean(j)+sqrt(2.0_dp)*scale*nodes))
         newv(j)=sum(w*(res%mean(j)+sqrt(2.0_dp)*scale*nodes-newm(j))**2)
         lse=mlog+log(den)+log(sqrt(2.0_dp)*scale); res%log_normalizer=res%log_normalizer+lse
      end do
      change=max(maxval(abs(newm-res%mean)),maxval(abs(newv-[(res%covariance(j,j),j=1,p)])))
      res%mean=newm; res%covariance=0.0_dp
      do j=1,p; res%covariance(j,j)=max(newv(j),1.0e-10_dp); end do
      res%iterations=it
      if(change<tol) then; res%converged=.true.; exit; end if
   end do
end subroutine componentwise_iterative_quadrature

subroutine adaptive_sparse_grid_quadrature(f,mean0,cov0,level,max_iter,tol,res)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: mean0(:),cov0(:,:)
   integer, intent(in) :: level,max_iter
   real(dp), intent(in) :: tol
   type(iterative_quad_result_t), intent(out) :: res
   integer :: p,it,info
   real(dp), allocatable :: l(:,:),mnew(:),cnew(:,:)
   real(dp) :: zscaled,mass,change,logjac,maxlog
   p=size(mean0)
   if(p>8) error stop 'adaptive_sparse_grid_quadrature: dimensions > 8 are unsupported'
   allocate(res%mean(p),res%covariance(p,p),l(p,p),mnew(p),cnew(p,p))
   res%mean=mean0; res%covariance=cov0; call make_positive_definite(res%covariance)
   do it=1,max_iter
      call chol_lower(res%covariance,l,info); if(info/=0) exit
      call sparse_integrals(f,res%mean,l,max(2,level),mass,mnew,cnew,maxlog)
      if(mass<=0.0_dp .or. mass/=mass) exit
      mnew=mnew/mass; cnew=cnew/mass-outer_product(mnew,mnew)
      cnew=cnew+1.0e-10_dp*identity_local(p); call make_positive_definite(cnew)
      change=max(maxval(abs(mnew-res%mean)),maxval(abs(cnew-res%covariance)))
      logjac=0.5_dp*real(p,dp)*log(2.0_dp)
      logjac=logjac+sum(log([(l(info,info),info=1,p)]))
      res%log_normalizer=log(abs(mass))+maxlog+logjac
      res%mean=mnew; res%covariance=cnew; res%iterations=it
      if(change<tol) then; res%converged=.true.; exit; end if
   end do
end subroutine adaptive_sparse_grid_quadrature

subroutine sparse_integrals(f,mu,l,q,mass,m1,m2,maxlog)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: mu(:),l(:,:)
   integer, intent(in) :: q
   real(dp), intent(out) :: mass,m1(:),m2(:,:),maxlog
   integer :: p
   integer, allocatable :: levels(:)
   real(dp) :: smass
   real(dp), allocatable :: sm1(:),sm2(:,:)
   p=size(mu); allocate(levels(p),sm1(p),sm2(p,p)); mass=0.0_dp; m1=0.0_dp; m2=0.0_dp
   maxlog=-huge(1.0_dp); call find_sparse_maxlog(f,mu,l,q,levels,1,maxlog)
   call accumulate_levels(f,mu,l,q,levels,1,maxlog,mass,m1,m2)
end subroutine sparse_integrals

recursive subroutine find_sparse_maxlog(f,mu,l,q,levels,pos,maxlog)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: mu(:),l(:,:)
   integer, intent(in) :: q,pos
   integer, intent(inout) :: levels(:)
   real(dp), intent(inout) :: maxlog
   integer :: lev,s,p,d
   p=size(levels); d=p
   if(pos>p) then
      s=sum(levels)
      if(s>=q-d+1 .and. s<=q) call tensor_find_max(f,mu,l,levels,maxlog)
      return
   end if
   do lev=1,q
      levels(pos)=lev
      if(sum(levels(1:pos))+(p-pos)<=q) call find_sparse_maxlog(f,mu,l,q,levels,pos+1,maxlog)
   end do
end subroutine find_sparse_maxlog

recursive subroutine accumulate_levels(f,mu,l,q,levels,pos,maxlog,mass,m1,m2)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: mu(:),l(:,:),maxlog
   integer, intent(in) :: q,pos
   integer, intent(inout) :: levels(:)
   real(dp), intent(inout) :: mass,m1(:),m2(:,:)
   integer :: lev,s,p,d,r,coef
   p=size(levels); d=p
   if(pos>p) then
      s=sum(levels)
      if(s<q-d+1 .or. s>q) return
      r=q-s; coef=(-1)**r*binomial_int(d-1,r)
      if(coef/=0) call tensor_accumulate(f,mu,l,levels,real(coef,dp),maxlog,mass,m1,m2)
      return
   end if
   do lev=1,q
      levels(pos)=lev
      if(sum(levels(1:pos))+(p-pos)<=q) call accumulate_levels(f,mu,l,q,levels,pos+1,maxlog,mass,m1,m2)
   end do
end subroutine accumulate_levels

subroutine tensor_find_max(f,mu,l,levels,maxlog)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: mu(:),l(:,:)
   integer, intent(in) :: levels(:)
   real(dp), intent(inout) :: maxlog
   integer :: p,d,total,i,tmp,idx,ord
   integer, allocatable :: orders(:)
   real(dp), allocatable :: nodes(:,:),weights(:,:),z(:),theta(:)
   real(dp) :: lv
   p=size(levels); allocate(orders(p)); orders=2*levels-1; total=product(orders)
   allocate(nodes(maxval(orders),p),weights(maxval(orders),p),z(p),theta(p)); nodes=0.0_dp; weights=0.0_dp
   do d=1,p; ord=orders(d); call gauss_hermite_rule(ord,nodes(1:ord,d),weights(1:ord,d)); end do
   do i=0,total-1
      tmp=i
      do d=1,p; idx=mod(tmp,orders(d))+1; tmp=tmp/orders(d); z(d)=nodes(idx,d); end do
      theta=mu+sqrt(2.0_dp)*matmul(l,z); lv=f(theta)+sum(z*z); if(lv>maxlog) maxlog=lv
   end do
end subroutine tensor_find_max

subroutine tensor_accumulate(f,mu,l,levels,coef,maxlog,mass,m1,m2)
   procedure(log_target_iface) :: f
   real(dp), intent(in) :: mu(:),l(:,:),coef,maxlog
   integer, intent(in) :: levels(:)
   real(dp), intent(inout) :: mass,m1(:),m2(:,:)
   integer :: p,d,total,i,tmp,idx,ord
   integer, allocatable :: orders(:)
   real(dp), allocatable :: nodes(:,:),weights(:,:),z(:),theta(:)
   real(dp) :: lv,wprod,w
   p=size(levels); allocate(orders(p)); orders=2*levels-1; total=product(orders)
   allocate(nodes(maxval(orders),p),weights(maxval(orders),p),z(p),theta(p)); nodes=0.0_dp; weights=0.0_dp
   do d=1,p; ord=orders(d); call gauss_hermite_rule(ord,nodes(1:ord,d),weights(1:ord,d)); end do
   do i=0,total-1
      tmp=i; wprod=coef
      do d=1,p
         idx=mod(tmp,orders(d))+1; tmp=tmp/orders(d); z(d)=nodes(idx,d); wprod=wprod*weights(idx,d)
      end do
      theta=mu+sqrt(2.0_dp)*matmul(l,z); lv=f(theta)+sum(z*z); w=wprod*exp(lv-maxlog)
      mass=mass+w; m1=m1+w*theta; m2=m2+w*outer_product(theta,theta)
   end do
end subroutine tensor_accumulate

pure integer function binomial_int(n,k) result(v)
   integer, intent(in) :: n,k
   integer :: i
   if(k<0 .or. k>n) then; v=0; return; end if
   v=1; do i=1,min(k,n-k); v=v*(n-i+1)/i; end do
end function binomial_int

pure function identity_local(n) result(a)
   integer, intent(in) :: n
   real(dp) :: a(n,n)
   integer :: i
   a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
end function identity_local

end module ld_quadrature_catalog
