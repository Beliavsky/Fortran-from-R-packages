module lavaan_browne
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : inverse_general, vech
   implicit none
   private

   type, public :: browne_test_result
      real(dp) :: statistic = huge(1.0_dp)
      real(dp) :: df = 0.0_dp
      integer :: status = 0
   end type browne_test_result

   public :: browne_residual_test, browne_residual_nt, browne_nt_gamma

contains

   subroutine browne_residual_test(residual, gamma, nobs, npar, result)
      real(dp), intent(in) :: residual(:), gamma(:, :)
      integer, intent(in) :: nobs, npar
      type(browne_test_result), intent(out) :: result
      real(dp), allocatable :: gi(:, :)
      integer :: info, q
      q = size(residual)
      if (nobs <= 0 .or. npar < 0 .or. npar >= q .or. any(shape(gamma) /= [q, q])) then
         result%status = -1
         return
      end if
      call inverse_general(gamma, gi, info)
      if (info /= 0) then
      result%status = info
      return
      end if
      result%statistic = real(nobs, dp) * dot_product(residual, matmul(gi, residual))
      result%df = real(q - npar, dp)
      result%status = 0
   end subroutine browne_residual_test

   subroutine browne_residual_nt(sample_cov, model_cov, nobs, npar, result)
      real(dp), intent(in) :: sample_cov(:, :), model_cov(:, :)
      integer, intent(in) :: nobs, npar
      type(browne_test_result), intent(out) :: result
      real(dp), allocatable :: gamma(:, :), residual(:)
      integer :: p
      p = size(sample_cov, 1)
      if (size(sample_cov, 2) /= p .or. any(shape(model_cov) /= [p, p])) then
         result%status = -1
         return
      end if
      call browne_nt_gamma(model_cov, gamma)
      residual = vech(sample_cov - model_cov)
      call browne_residual_test(residual, gamma, nobs, npar, result)
   end subroutine browne_residual_nt

   subroutine browne_nt_gamma(sigma, gamma)
      real(dp), intent(in) :: sigma(:, :)
      real(dp), allocatable, intent(out) :: gamma(:, :)
      integer :: p, q, i, j, k, l, a, b
      p = size(sigma, 1)
      q = p * (p + 1) / 2
      allocate(gamma(q, q))
      gamma = 0.0_dp
      a = 0
      do j = 1, p
         do i = j, p
            a = a + 1
            b = 0
            do l = 1, p
               do k = l, p
                  b = b + 1
                  gamma(a, b) = sigma(i, k) * sigma(j, l) + sigma(i, l) * sigma(j, k)
               end do
            end do
         end do
      end do
   end subroutine browne_nt_gamma

end module lavaan_browne
