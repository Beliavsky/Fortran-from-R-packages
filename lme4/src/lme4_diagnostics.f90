module lme4_diagnostics
   use lme4_kinds, only : dp
   use lme4_types, only : covariance_block_t, lmm_result_t, glmm_result_t, &
      family_binomial, family_poisson, family_gamma, family_inverse_gaussian, &
      family_negative_binomial
   use lme4_covariance, only : covariance_pca
   implicit none
   private
   public :: is_singular, re_pca, standardized_residuals_lmm, pearson_residuals_glmm

   interface is_singular
      module procedure is_singular_lmm
      module procedure is_singular_glmm
   end interface

contains

   logical function is_singular_lmm(result,tolerance) result(singular)
      type(lmm_result_t), intent(in) :: result
      real(dp), intent(in), optional :: tolerance
      singular = blocks_singular(result%varcorr,tolerance)
   end function is_singular_lmm

   logical function is_singular_glmm(result,tolerance) result(singular)
      type(glmm_result_t), intent(in) :: result
      real(dp), intent(in), optional :: tolerance
      singular = blocks_singular(result%varcorr,tolerance)
   end function is_singular_glmm

   logical function blocks_singular(blocks,tolerance) result(singular)
      type(covariance_block_t), intent(in) :: blocks(:)
      real(dp), intent(in), optional :: tolerance
      real(dp) :: tol
      real(dp), allocatable :: sd(:), rotation(:,:)
      integer :: k, info
      tol = 1.0e-4_dp
      if (present(tolerance)) tol = tolerance
      singular = .false.
      do k = 1, size(blocks)
         call covariance_pca(blocks(k)%covariance,sd,rotation,info)
         if (info /= 0 .or. minval(sd) < tol*max(1.0_dp,maxval(sd))) then
            singular = .true.
            return
         end if
      end do
   end function blocks_singular

   subroutine re_pca(block, standard_deviations, rotation, info)
      type(covariance_block_t), intent(in) :: block
      real(dp), allocatable, intent(out) :: standard_deviations(:), rotation(:,:)
      integer, intent(out) :: info
      call covariance_pca(block%covariance,standard_deviations,rotation,info)
   end subroutine re_pca

   subroutine standardized_residuals_lmm(result,residuals)
      type(lmm_result_t), intent(in) :: result
      real(dp), allocatable, intent(out) :: residuals(:)
      if (result%sigma > 0.0_dp) then
         residuals = result%residuals/result%sigma
      else
         residuals = result%residuals
      end if
   end subroutine standardized_residuals_lmm

   subroutine pearson_residuals_glmm(result,y,residuals)
      type(glmm_result_t), intent(in) :: result
      real(dp), intent(in) :: y(:)
      real(dp), allocatable, intent(out) :: residuals(:)
      real(dp), allocatable :: variance(:)
      allocate(variance(size(y)))
      select case (result%family)
      case (family_binomial)
         variance = result%fitted*(1.0_dp-result%fitted)
      case (family_poisson)
         variance = result%fitted
      case (family_gamma)
         variance = result%dispersion*result%fitted**2
      case (family_inverse_gaussian)
         variance = result%dispersion*result%fitted**3
      case (family_negative_binomial)
         variance = result%fitted+result%fitted**2/result%dispersion
      case default
         variance = 1.0_dp
      end select
      residuals = (y-result%fitted)/sqrt(max(variance,1.0e-12_dp))
   end subroutine pearson_residuals_glmm

end module lme4_diagnostics
