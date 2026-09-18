module gamm4_covariance
   use gamm4_kinds, only : dp
   use gamm4_types, only : gamm4_vb_result_t
   use r_linalg, only : spd_inverse_logdet, cholesky_factor
   implicit none
   private
   public :: gamm4_get_vb

contains

   subroutine gamm4_get_vb(v, z, phi, scale, xf, xfp, sp, b, result)
      real(dp), intent(in) :: v(:) !! Sampling-variance diagonal before adding non-smooth random-effect covariance.
      real(dp), intent(in) :: z(:, :) !! Expanded non-smooth random-effect design matrix with one row per observation.
      real(dp), intent(in) :: phi(:, :) !! Relative covariance matrix for non-smooth random coefficients.
      real(dp), intent(in) :: scale !! Response scale multiplying the relative random-effect covariance.
      real(dp), intent(in) :: xf(:, :) !! GAM design matrix in the original smooth parameterization.
      real(dp), intent(in) :: xfp(:, :) !! GAM design matrix in the fitting/reparameterized smooth parameterization.
      real(dp), intent(in) :: sp(:) !! Smoothing parameters aligned with columns of `xfp`; zero means unpenalized.
      real(dp), intent(in) :: b(:, :) !! Linear map from fitting-parameter coefficients to original GAM coefficients.
      type(gamm4_vb_result_t), intent(out) :: result !! Marginal coefficient covariance, information, and Cholesky factor.
      real(dp), allocatable :: vinv(:, :), vmat(:, :), xvxs(:, :), xvxs_inv(:, :)
      real(dp) :: ignored_logdet
      integer :: i, info, n, p

      result%status = 0
      n = size(v)
      p = size(xf, 2)
      if (scale <= 0.0_dp .or. size(xf, 1) /= n .or. size(xfp, 1) /= n .or. &
          size(xfp, 2) /= p .or. size(sp) /= p .or. size(b, 1) /= p .or. size(b, 2) /= p) then
         result%status = 1
         return
      end if
      if (size(z, 1) /= n .or. size(phi, 1) /= size(z, 2) .or. size(phi, 2) /= size(z, 2)) then
         result%status = 2
         return
      end if
      if (any(v <= 0.0_dp) .or. any(sp < 0.0_dp)) then
         result%status = 3
         return
      end if

      allocate(vmat(n, n))
      vmat = 0.0_dp
      do i = 1, n
         vmat(i, i) = v(i)
      end do
      if (size(z, 2) > 0) vmat = vmat + scale * matmul(matmul(z, phi), transpose(z))
      call spd_inverse_logdet(vmat, vinv, ignored_logdet, info)
      if (info /= 0) then
         result%status = 4
         return
      end if
      result%xvx = matmul(transpose(xf), matmul(vinv, xf))
      xvxs = matmul(transpose(xfp), matmul(vinv, xfp))
      do i = 1, p
         xvxs(i, i) = xvxs(i, i) + sp(i) / scale
      end do
      call spd_inverse_logdet(xvxs, xvxs_inv, ignored_logdet, info)
      if (info /= 0) then
         result%status = 5
         return
      end if
      result%vb = matmul(matmul(b, xvxs_inv), transpose(b))
      call cholesky_factor(result%xvx, result%r, info)
      if (info /= 0) then
         if (allocated(result%r)) deallocate(result%r)
         allocate(result%r(0, 0))
      end if
   end subroutine gamm4_get_vb

end module gamm4_covariance
