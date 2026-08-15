! Distribution-comparison helpers corresponding to fitDist/chooseDist computational work.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_selection
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar
   use gamlss_types
   use gamlss_core, only : fit_gamlss_model
   implicit none
   private
   public :: gamlss_family_comparison_t, compare_families, best_family

   type, public :: gamlss_family_comparison_t
      integer, allocatable :: family(:),status(:)
      real(dp), allocatable :: global_deviance(:),df(:),gaic(:),sbc(:)
   end type gamlss_family_comparison_t

contains

   subroutine compare_families(y,families,comparison,k,weights,control)
      real(dp),intent(in)::y(:)
      integer,intent(in)::families(:)
      type(gamlss_family_comparison_t),intent(out)::comparison
      real(dp),intent(in),optional::k,weights(:)
      type(gamlss_control_t),intent(in),optional::control
      real(dp),allocatable::x(:,:)
      type(gamlss_result_t)::fit
      real(dp)::kk
      integer::i,np
      allocate(x(size(y),1));x=1.0_dp;kk=2.0_dp;if(present(k))kk=k
      allocate(comparison%family(size(families)),comparison%status(size(families)))
      allocate(comparison%global_deviance(size(families)),comparison%df(size(families)))
      allocate(comparison%gaic(size(families)),comparison%sbc(size(families)))
      do i=1,size(families)
         np=family_npar(families(i));comparison%family(i)=families(i)
         select case(np)
         case(1)
            call fit_gamlss_model(y,x,families(i),fit,weights=weights,control=control)
         case(2)
            call fit_gamlss_model(y,x,families(i),fit,x_sigma=x,weights=weights,control=control)
         case(3)
            call fit_gamlss_model(y,x,families(i),fit,x_sigma=x,x_nu=x,weights=weights,control=control)
         case(4)
            call fit_gamlss_model(y,x,families(i),fit,x_sigma=x,x_nu=x,x_tau=x,weights=weights,control=control)
         case default
            fit%status=99
         end select
         comparison%status(i)=fit%status
         if(fit%status==0)then
            comparison%global_deviance(i)=fit%global_deviance;comparison%df(i)=fit%df_fit
            comparison%gaic(i)=fit%global_deviance+kk*fit%df_fit;comparison%sbc(i)=fit%sbc
         else
            comparison%global_deviance(i)=huge(1.0_dp);comparison%df(i)=0.0_dp
            comparison%gaic(i)=huge(1.0_dp);comparison%sbc(i)=huge(1.0_dp)
         end if
      end do
   end subroutine compare_families

   integer function best_family(comparison,use_sbc) result(family)
      type(gamlss_family_comparison_t),intent(in)::comparison
      logical,intent(in),optional::use_sbc
      logical::bic
      integer::idx
      bic=.false.;if(present(use_sbc))bic=use_sbc
      if(size(comparison%family)==0)then;family=0;return;end if
      if(bic)then;idx=minloc(comparison%sbc,dim=1);else;idx=minloc(comparison%gaic,dim=1);end if
      family=comparison%family(idx)
   end function best_family

end module gamlss_selection
