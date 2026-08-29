! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_miph
   use r_compat, only: dp
   use matrixdist_iph, only: transform_time_jac
   use matrixdist_multivariate, only: bivph_density, bivph_tail
   use matrixdist_ph, only: ph_density, ph_cdf
   implicit none
   private
   public :: biviph_density, biviph_tail
   public :: miph_density_point, miph_cdf_point
contains
   function biviph_density(x, alpha, s11, s12, s22, kinds, beta) result(f)
      real(dp), intent(in) :: x(2), alpha(:), s11(:,:), s12(:,:), s22(:,:), beta(:,:)
      character(len=*), intent(in) :: kinds(2)
      real(dp) :: f, t1, t2, j1, j2
      logical :: ok1, ok2
      call transform_time_jac(x(1), kinds(1), beta(:,1), t1, j1, ok1)
      call transform_time_jac(x(2), kinds(2), beta(:,2), t2, j2, ok2)
      if (.not.(ok1.and.ok2)) then
         f=0.0_dp
      else
         f=bivph_density(t1,t2,alpha,s11,s12,s22)*j1*j2
      end if
   end function biviph_density

   function biviph_tail(x, alpha, s11, s12, s22, kinds, beta) result(q)
      real(dp), intent(in) :: x(2), alpha(:), s11(:,:), s12(:,:), s22(:,:), beta(:,:)
      character(len=*), intent(in) :: kinds(2)
      real(dp) :: q, t1, t2, j1, j2
      logical :: ok1, ok2
      call transform_time_jac(x(1), kinds(1), beta(:,1), t1, j1, ok1)
      call transform_time_jac(x(2), kinds(2), beta(:,2), t2, j2, ok2)
      if (.not.(ok1.and.ok2)) then
         q=0.0_dp
      else
         q=bivph_tail(t1,t2,alpha,s11,s12,s22)
      end if
   end function biviph_tail

   function miph_density_point(y,alpha,s,kinds,beta,delta) result(v)
      real(dp),intent(in)::y(:),alpha(:),s(:,:,:),beta(:,:)
      character(len=*),intent(in)::kinds(:)
      logical,intent(in),optional::delta(:)
      real(dp)::v,prodv,t,jac,f
      real(dp),allocatable::unit(:)
      integer::p,d,j,i
      logical::ok,unc
      p=size(alpha)
      d=size(s,3)
      if(size(y)/=d .or. size(kinds)/=d .or. size(beta,2)/=d) error stop 'miph_density_point: dimension mismatch'
      allocate(unit(p))
      v=0.0_dp
      do j=1,p
         unit=0.0_dp
         unit(j)=1.0_dp
         prodv=1.0_dp
         do i=1,d
            call transform_time_jac(y(i),kinds(i),beta(:,i),t,jac,ok)
            if(.not.ok) then
            prodv=0.0_dp
            exit
            end if
            unc=.true.
            if(present(delta))unc=delta(i)
            if(unc) then
               f=ph_density(t,unit,s(:,:,i))*jac
            else
               f=1.0_dp-ph_cdf(t,unit,s(:,:,i))
            end if
            prodv=prodv*f
         end do
         v=v+alpha(j)*prodv
      end do
   end function miph_density_point

   function miph_cdf_point(y,alpha,s,kinds,beta,lower_tail) result(v)
      real(dp),intent(in)::y(:),alpha(:),s(:,:,:),beta(:,:)
      character(len=*),intent(in)::kinds(:)
      logical,intent(in),optional::lower_tail
      real(dp)::v,prodv,t,jac
      real(dp),allocatable::unit(:)
      integer::p,d,j,i
      logical::ok,lower
      lower=.true.
      if(present(lower_tail))lower=lower_tail
      p=size(alpha)
      d=size(s,3)
      if(size(y)/=d .or. size(kinds)/=d .or. size(beta,2)/=d) error stop 'miph_cdf_point: dimension mismatch'
      allocate(unit(p))
      v=0.0_dp
      do j=1,p
         unit=0.0_dp
         unit(j)=1.0_dp
         prodv=1.0_dp
         do i=1,d
            call transform_time_jac(y(i),kinds(i),beta(:,i),t,jac,ok)
            if(.not.ok) then
               if(lower) then
               prodv=0.0_dp
               else
               prodv=prodv
               end if
               cycle
            end if
            prodv=prodv*ph_cdf(t,unit,s(:,:,i),lower)
         end do
         v=v+alpha(j)*prodv
      end do
   end function miph_cdf_point
end module matrixdist_miph
