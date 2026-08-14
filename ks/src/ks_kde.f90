! SPDX-License-Identifier: GPL-2.0-only
module ks_kde
   use ks_kinds, only: dp, pi, log2pi
   use ks_linalg, only: covariance_matrix, spd_inverse, spd_logdet, spd_cholesky
   use ks_normal, only: normal_cdf, normal_pdf, mvn_derivative_tensor
   use ks_rng, only: rng_state, rng_uniform, rng_normal
   implicit none
   private
   public :: kde_model, fit_kde, kde_pdf, kde_logpdf, kdde_eval
   public :: kde_cdf_1d, kde_quantile_1d, kde_random, kfe_tensor
   public :: contour_levels_grid, contour_probability_grid, truncate_grid_density
   public :: positive_kde_pdf_1d, unit_interval_kde_pdf_1d

   type :: kde_model
      integer :: n=0,d=0
      real(dp), allocatable :: x(:,:),w(:),H(:,:),Hinv(:,:),cholH(:,:)
      real(dp) :: log_norm=0.0_dp
   end type kde_model
contains
   subroutine default_hns(x,h,deriv_order)
      real(dp),intent(in)::x(:,:)
      real(dp),intent(out)::h(size(x,2),size(x,2))
      integer,intent(in),optional::deriv_order
      real(dp)::s(size(x,2),size(x,2)),fac
      integer::r,n,d
      r=0;if(present(deriv_order))r=deriv_order
      n=size(x,1);d=size(x,2);call covariance_matrix(x,s)
      fac=(4.0_dp/(real(n,dp)*real(d+2*r+2,dp)))**(2.0_dp/real(d+2*r+4,dp))
      h=fac*s
   end subroutine default_hns

   subroutine fit_kde(x,model,H,weights,info)
      real(dp),intent(in)::x(:,:)
      type(kde_model),intent(out)::model
      real(dp),intent(in),optional::H(:,:),weights(:)
      integer,intent(out),optional::info
      integer::ierr
      real(dp)::sw,ld
      model%n=size(x,1);model%d=size(x,2)
      allocate(model%x(model%n,model%d),model%w(model%n),model%H(model%d,model%d), &
               model%Hinv(model%d,model%d),model%cholH(model%d,model%d))
      model%x=x
      if(present(weights)) then
         model%w=weights
         sw=sum(model%w)
         if(sw>0.0_dp) model%w=model%w*real(model%n,dp)/sw
      else
         model%w=1.0_dp
      end if
      if(present(H)) then; model%H=H; else; call default_hns(x,model%H); end if
      call spd_inverse(model%H,model%Hinv,ierr)
      if(ierr==0) call spd_cholesky(model%H,model%cholH,ierr)
      if(ierr==0) then
         ld=spd_logdet(model%H,ierr)
         model%log_norm=-0.5_dp*(real(model%d,dp)*log2pi+ld)
      else
         model%log_norm=-huge(1.0_dp)
      end if
      if(present(info))info=ierr
   end subroutine fit_kde

   subroutine kde_pdf(model,points,values)
      type(kde_model),intent(in)::model
      real(dp),intent(in)::points(:,:)
      real(dp),intent(out)::values(size(points,1))
      integer::i,j
      real(dp)::dx(model%d),q
      values=0.0_dp
      do i=1,size(points,1)
         do j=1,model%n
            dx=points(i,:)-model%x(j,:)
            q=dot_product(dx,matmul(model%Hinv,dx))
            if(q<1500.0_dp) values(i)=values(i)+model%w(j)*exp(model%log_norm-0.5_dp*q)
         end do
         values(i)=values(i)/real(model%n,dp)
      end do
   end subroutine kde_pdf

   subroutine kde_logpdf(model,points,values)
      type(kde_model),intent(in)::model
      real(dp),intent(in)::points(:,:)
      real(dp),intent(out)::values(size(points,1))
      real(dp)::dens(size(points,1))
      call kde_pdf(model,points,dens)
      values=log(max(dens,tiny(1.0_dp)))
   end subroutine kde_logpdf

   subroutine kdde_eval(model,points,deriv_order,deriv)
      type(kde_model),intent(in)::model
      real(dp),intent(in)::points(:,:)
      integer,intent(in)::deriv_order
      real(dp),allocatable,intent(out)::deriv(:,:)
      integer::i,j,m
      real(dp),allocatable::tmp(:)
      m=model%d**deriv_order;allocate(deriv(size(points,1),m));deriv=0.0_dp
      do i=1,size(points,1)
         do j=1,model%n
            call mvn_derivative_tensor(points(i,:),model%x(j,:),model%H,deriv_order,tmp)
            deriv(i,:)=deriv(i,:)+model%w(j)*tmp
         end do
         deriv(i,:)=deriv(i,:)/real(model%n,dp)
      end do
   end subroutine kdde_eval

   function kde_cdf_1d(model,x) result(v)
      type(kde_model),intent(in)::model
      real(dp),intent(in)::x
      real(dp)::v,h
      integer::j
      if(model%d/=1) then;v=0.0_dp;return;end if
      h=sqrt(model%H(1,1));v=0.0_dp
      do j=1,model%n
         v=v+model%w(j)*normal_cdf(x,model%x(j,1),h)
      end do
      v=v/real(model%n,dp)
   end function kde_cdf_1d

   function kde_quantile_1d(model,p,tol) result(q)
      type(kde_model),intent(in)::model
      real(dp),intent(in)::p
      real(dp),intent(in),optional::tol
      real(dp)::q,lo,hi,mid,eps,h
      integer::it
      eps=1e-10_dp;if(present(tol))eps=tol
      if(model%d/=1)then;q=0.0_dp;return;end if
      h=sqrt(model%H(1,1));lo=minval(model%x(:,1))-10*h;hi=maxval(model%x(:,1))+10*h
      do it=1,200
         mid=0.5_dp*(lo+hi)
         if(kde_cdf_1d(model,mid)<p)then;lo=mid;else;hi=mid;end if
         if(abs(hi-lo)<=eps*max(1.0_dp,abs(mid)))exit
      end do
      q=0.5_dp*(lo+hi)
   end function kde_quantile_1d

   subroutine kde_random(model,rng,sample)
      type(kde_model),intent(in)::model
      type(rng_state),intent(inout)::rng
      real(dp),intent(out)::sample(:,:)
      integer::i,j,k
      real(dp)::u,c,z(model%d)
      do i=1,size(sample,1)
         u=rng_uniform(rng)*real(model%n,dp);c=0.0_dp;k=model%n
         do j=1,model%n
            c=c+model%w(j)
            if(u<=c)then;k=j;exit;end if
         end do
         do j=1,model%d;z(j)=rng_normal(rng);end do
         sample(i,:)=model%x(k,:)+matmul(model%cholH,z)
      end do
   end subroutine kde_random

   subroutine kfe_tensor(x,G,deriv_order,inc,estimate,weights)
      real(dp),intent(in)::x(:,:),G(:,:)
      integer,intent(in)::deriv_order
      logical,intent(in)::inc
      real(dp),allocatable,intent(out)::estimate(:)
      real(dp),intent(in),optional::weights(:)
      integer::n,d,m,i,j,npair
      real(dp)::wi,wj,denom
      real(dp),allocatable::tmp(:)
      real(dp)::zero(size(x,2))
      n=size(x,1);d=size(x,2);m=d**deriv_order;allocate(estimate(m));estimate=0.0_dp;zero=0.0_dp
      denom=0.0_dp;npair=0
      do i=1,n
         wi=1.0_dp;if(present(weights))wi=weights(i)
         do j=1,n
            if(.not.inc .and. i==j)cycle
            wj=1.0_dp;if(present(weights))wj=weights(j)
            call mvn_derivative_tensor(x(i,:)-x(j,:),zero,G,deriv_order,tmp)
            estimate=estimate+wi*wj*tmp;denom=denom+wi*wj;npair=npair+1
         end do
      end do
      if(denom>0.0_dp)estimate=estimate/denom
   end subroutine kfe_tensor

   subroutine sort_desc_with_index(a,idx)
      real(dp),intent(in)::a(:)
      integer,intent(out)::idx(size(a))
      integer::i,j,key
      do i=1,size(a);idx(i)=i;end do
      do i=2,size(a)
         key=idx(i);j=i-1
         do while(j>=1)
            if(a(idx(j))>=a(key))exit
            idx(j+1)=idx(j);j=j-1
         end do
         idx(j+1)=key
      end do
   end subroutine sort_desc_with_index

   subroutine contour_levels_grid(density,cell_volume,probs,levels)
      real(dp),intent(in)::density(:),cell_volume,probs(:)
      real(dp),intent(out)::levels(size(probs))
      integer,allocatable::idx(:)
      real(dp)::total,cum,target
      integer::p,k
      allocate(idx(size(density)));call sort_desc_with_index(density,idx)
      total=sum(max(density,0.0_dp))*cell_volume
      do p=1,size(probs)
         target=max(0.0_dp,min(1.0_dp,probs(p)))*total;cum=0.0_dp;levels(p)=0.0_dp
         do k=1,size(idx)
            cum=cum+max(0.0_dp,density(idx(k)))*cell_volume
            if(cum>=target)then;levels(p)=density(idx(k));exit;end if
         end do
      end do
   end subroutine contour_levels_grid

   function contour_probability_grid(density,cell_volume,level) result(prob)
      real(dp),intent(in)::density(:),cell_volume,level
      real(dp)::prob,total
      total=sum(max(density,0.0_dp))*cell_volume
      if(total<=0.0_dp)then;prob=0.0_dp;else
         prob=sum(merge(max(density,0.0_dp),0.0_dp,density>=level))*cell_volume/total
      end if
   end function contour_probability_grid

   subroutine truncate_grid_density(density,inside,cell_volume,out)
      real(dp),intent(in)::density(:),cell_volume
      logical,intent(in)::inside(:)
      real(dp),intent(out)::out(size(density))
      real(dp)::s0,s1
      s0=sum(max(density,0.0_dp))*cell_volume
      out=merge(density,0.0_dp,inside)
      s1=sum(max(out,0.0_dp))*cell_volume
      if(s1>0.0_dp)out=out*s0/s1
   end subroutine truncate_grid_density

   function positive_kde_pdf_1d(data,h,x,adj,weights) result(v)
      real(dp),intent(in)::data(:),h,x
      real(dp),intent(in),optional::adj,weights(:)
      real(dp)::v,a,y,yi,w,sw
      integer::i
      a=abs(minval(data));if(present(adj))a=adj
      if(x+a<=0.0_dp)then;v=0.0_dp;return;end if
      y=log(x+a);v=0.0_dp;sw=0.0_dp
      do i=1,size(data)
         if(data(i)+a<=0.0_dp)cycle
         yi=log(data(i)+a);w=1.0_dp;if(present(weights))w=weights(i)
         v=v+w*normal_pdf(y,yi,h);sw=sw+w
      end do
      if(sw>0.0_dp)v=v/sw/(x+a)
   end function positive_kde_pdf_1d

   function unit_interval_kde_pdf_1d(data,h,x,weights) result(v)
      real(dp),intent(in)::data(:),h,x
      real(dp),intent(in),optional::weights(:)
      real(dp)::v,y,yi,w,sw,jac
      integer::i
      if(x<=0.0_dp .or. x>=1.0_dp)then;v=0.0_dp;return;end if
      y=log(x/(1.0_dp-x));jac=1.0_dp/(x*(1.0_dp-x));v=0.0_dp;sw=0.0_dp
      do i=1,size(data)
         if(data(i)<=0.0_dp .or. data(i)>=1.0_dp)cycle
         yi=log(data(i)/(1.0_dp-data(i)));w=1.0_dp;if(present(weights))w=weights(i)
         v=v+w*normal_pdf(y,yi,h);sw=sw+w
      end do
      if(sw>0.0_dp)v=v/sw*jac
   end function unit_interval_kde_pdf_1d
end module ks_kde
