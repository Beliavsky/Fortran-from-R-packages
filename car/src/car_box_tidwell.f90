module car_box_tidwell
   use car_kinds, only: dp
   use car_utils, only: f_survival, least_squares_residuals, matrix_inverse, student_t_two_sided_p
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   implicit none
   private

   public :: box_tidwell

contains

   subroutine box_tidwell(y, x1, powers, score_t, p_values, iterations, status, x2, max_iter, tol, &
                          omnibus_f, omnibus_p, converged)
      real(dp), intent(in) :: y(:) !! Response vector with one observation per row of x1 and optional x2.
      real(dp), intent(in) :: x1(:, :) !! Strictly positive predictor columns whose powers are to be estimated.
      real(dp), allocatable, intent(out) :: powers(:) !! Estimated Box-Tidwell powers, one for each x1 column.
      real(dp), allocatable, intent(out) :: score_t(:) !! Initial score t statistics for testing unit powers.
      real(dp), allocatable, intent(out) :: p_values(:) !! Two-sided Student-t probabilities for score_t.
      integer, intent(out) :: iterations !! Number of power-update iterations performed.
      integer, intent(out) :: status !! Zero on a usable fit; nonzero for invalid data, singular fits, or numerical failure.
      real(dp), intent(in), optional :: x2(:, :) !! Optional additional predictor columns that are not transformed.
      integer, intent(in), optional :: max_iter !! Maximum requested iterations; default 25, with R-compatible one-final-update semantics.
      real(dp), intent(in), optional :: tol !! Positive relative convergence tolerance; default 0.001.
      real(dp), intent(out), optional :: omnibus_f !! Optional initial omnibus F statistic for all transformation terms.
      real(dp), intent(out), optional :: omnibus_p !! Optional upper-tail F probability corresponding to omnibus_f.
      logical, intent(out), optional :: converged !! True when the relative power update met the tolerance before the iteration limit.

      real(dp), allocatable :: beta1(:)
      real(dp), allocatable :: beta2(:)
      real(dp), allocatable :: design1(:, :)
      real(dp), allocatable :: design2(:, :)
      real(dp), allocatable :: hcov(:, :)
      real(dp), allocatable :: hcov_inv(:, :)
      real(dp), allocatable :: inv_xtx2(:, :)
      real(dp), allocatable :: last_powers(:)
      real(dp), allocatable :: residual1(:)
      real(dp), allocatable :: residual2(:)
      real(dp), allocatable :: transformed(:, :)
      real(dp), allocatable :: xlogx(:, :)
      real(dp), allocatable :: xtx2(:, :)
      real(dp) :: f_stat
      real(dp) :: rel_change
      real(dp) :: sigma2
      real(dp) :: tolerance
      integer :: df2
      integer :: fit_status
      integer :: i
      integer :: inv_status
      integer :: k1
      integer :: k2
      integer :: limit
      integer :: n
      integer :: p1
      integer :: p2
      logical :: did_converge

      n = size(y)
      k1 = size(x1, 2)
      k2 = 0
      if (present(x2)) k2 = size(x2, 2)
      p1 = 1 + k1 + k2
      p2 = 1 + 2 * k1 + k2

      allocate(powers(k1), score_t(k1), p_values(k1))
      powers = ieee_value(0.0_dp, ieee_quiet_nan)
      score_t = ieee_value(0.0_dp, ieee_quiet_nan)
      p_values = ieee_value(0.0_dp, ieee_quiet_nan)
      iterations = 0
      status = 0
      did_converge = .false.
      if (present(omnibus_f)) omnibus_f = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(omnibus_p)) omnibus_p = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(converged)) converged = .false.

      limit = 25
      if (present(max_iter)) limit = max_iter
      tolerance = 0.001_dp
      if (present(tol)) tolerance = tol

      if (k1 < 1 .or. n /= size(x1, 1) .or. n <= p2 .or. limit < 0 .or. tolerance <= 0.0_dp) then
         status = 1
         return
      end if
      if (present(x2)) then
         if (size(x2, 1) /= n) then
            status = 1
            return
         end if
      end if
      if (any(x1 <= 0.0_dp) .or. any(ieee_is_nan(x1)) .or. any(.not. ieee_is_finite(x1)) .or. &
          any(ieee_is_nan(y)) .or. any(.not. ieee_is_finite(y))) then
         status = 2
         return
      end if
      if (present(x2)) then
         if (any(ieee_is_nan(x2)) .or. any(.not. ieee_is_finite(x2))) then
            status = 2
            return
         end if
      end if

      allocate(design1(n, p1), design2(n, p2), beta1(p1), beta2(p2), residual1(n), residual2(n))
      allocate(transformed(n, k1), xlogx(n, k1), last_powers(k1))

      transformed = x1
      xlogx = x1 * log(x1)
      call fill_designs(transformed, xlogx, design1, design2)
      call least_squares_residuals(design1, y, beta1, residual1, fit_status)
      if (fit_status /= 0) then
         status = 3
         return
      end if
      call least_squares_residuals(design2, y, beta2, residual2, fit_status)
      if (fit_status /= 0) then
         status = 3
         return
      end if

      df2 = n - p2
      sigma2 = sum(residual2 * residual2) / real(df2, dp)
      allocate(xtx2(p2, p2), inv_xtx2(p2, p2))
      xtx2 = matmul(transpose(design2), design2)
      call matrix_inverse(xtx2, inv_xtx2, inv_status)
      if (inv_status /= 0 .or. sigma2 < 0.0_dp .or. .not. ieee_is_finite(sigma2)) then
         status = 3
         return
      end if

      do i = 1, k1
         if (sigma2 * inv_xtx2(1 + i, 1 + i) <= 0.0_dp) then
            status = 4
            return
         end if
         score_t(i) = beta2(1 + i) / sqrt(sigma2 * inv_xtx2(1 + i, 1 + i))
         p_values(i) = student_t_two_sided_p(score_t(i), df2)
         if (abs(beta1(1 + i)) <= tiny(1.0_dp)) then
            status = 4
            return
         end if
         powers(i) = 1.0_dp + beta2(1 + i) / beta1(1 + i)
      end do

      if (present(omnibus_f) .or. present(omnibus_p)) then
         allocate(hcov(k1, k1), hcov_inv(k1, k1))
         hcov = sigma2 * inv_xtx2(2:1 + k1, 2:1 + k1)
         call matrix_inverse(hcov, hcov_inv, inv_status)
         if (inv_status /= 0) then
            status = 4
            return
         end if
         f_stat = dot_product(beta2(2:1 + k1), matmul(hcov_inv, beta2(2:1 + k1))) / real(k1, dp)
         if (present(omnibus_f)) omnibus_f = f_stat
         if (present(omnibus_p)) omnibus_p = f_survival(f_stat, k1, df2)
      end if

      last_powers = 1.0_dp
      rel_change = maxval(abs((powers - last_powers) / (powers + tolerance)))
      do while (rel_change > tolerance .and. iterations <= limit)
         iterations = iterations + 1
         do i = 1, k1
            transformed(:, i) = x1(:, i)**powers(i)
         end do
         if (any(.not. ieee_is_finite(transformed)) .or. any(transformed <= 0.0_dp)) then
            status = 5
            return
         end if
         xlogx = transformed * log(transformed)
         call fill_designs(transformed, xlogx, design1, design2)
         call least_squares_residuals(design1, y, beta1, residual1, fit_status)
         if (fit_status /= 0) then
            status = 3
            return
         end if
         call least_squares_residuals(design2, y, beta2, residual2, fit_status)
         if (fit_status /= 0) then
            status = 3
            return
         end if

         last_powers = powers
         do i = 1, k1
            if (abs(beta1(1 + i)) <= tiny(1.0_dp)) then
               status = 4
               return
            end if
            powers(i) = powers(i) * (1.0_dp + beta2(1 + i) / beta1(1 + i))
         end do
         if (any(.not. ieee_is_finite(powers))) then
            status = 5
            return
         end if
         rel_change = maxval(abs((powers - last_powers) / (powers + tolerance)))
      end do

      did_converge = rel_change <= tolerance
      if (present(converged)) converged = did_converge

   contains

      pure subroutine fill_designs(transformed_local, xlogx_local, design1_local, design2_local)
         real(dp), intent(in) :: transformed_local(:, :) !! Current transformed x1 matrix used in both regression designs.
         real(dp), intent(in) :: xlogx_local(:, :) !! Current transformed predictors multiplied by their logarithms.
         real(dp), intent(out) :: design1_local(:, :) !! Regression design containing intercept, transformed x1, and optional x2.
         real(dp), intent(out) :: design2_local(:, :) !! Augmented design containing transformation score columns before design1 terms.
         integer :: first_x2

         design1_local(:, 1) = 1.0_dp
         design1_local(:, 2:1 + k1) = transformed_local
         design2_local(:, 1) = 1.0_dp
         design2_local(:, 2:1 + k1) = xlogx_local
         design2_local(:, 2 + k1:1 + 2 * k1) = transformed_local
         if (k2 > 0) then
            first_x2 = 2 + k1
            design1_local(:, first_x2:first_x2 + k2 - 1) = x2
            first_x2 = 2 + 2 * k1
            design2_local(:, first_x2:first_x2 + k2 - 1) = x2
         end if
      end subroutine fill_designs

   end subroutine box_tidwell

end module car_box_tidwell
