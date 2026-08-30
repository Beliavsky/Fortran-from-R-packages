! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
module vars_causality
   use r_kinds, only : dp
   use r_distributions, only : r_pchisq, r_pf
   use r_linalg, only : inverse_matrix
   use vars_types
   use vars_utils, only : duplication_matrix, kronecker_product, vech_lower
   implicit none
   private

   public :: granger_causality, instantaneous_causality

contains

   subroutine granger_causality(model, cause, result, info)
      type(var_model), intent(in) :: model
      logical, intent(in) :: cause(:)
      type(vars_test_result), intent(out) :: result
      integer, intent(out) :: info
      real(dp), allocatable :: xtxinv(:, :), sigma(:, :), values(:), covariance(:, :), invcov(:, :)
      integer, allocatable :: response_idx(:), reg_idx(:)
      integer :: ncause, nother, nrest, a, b, eq, lag, cause_var, idx, invinfo
      real(dp) :: denom

      if (size(cause) /= model%k) then
         info = vars_invalid_argument
         return
      end if
      ncause = count(cause)
      nother = model%k - ncause
      if (ncause < 1 .or. nother < 1) then
         info = vars_invalid_argument
         return
      end if
      if (.not. all(model%active)) then
         info = vars_invalid_argument
         return
      end if
      call inverse_matrix(matmul(transpose(model%x), model%x), xtxinv, info)
      if (info /= 0) then
         info = vars_singular
         return
      end if
      denom = real(model%nobs - model%nreg, dp)
      allocate(sigma(model%k, model%k))
      sigma = matmul(transpose(model%resid), model%resid) / denom
      nrest = model%p * ncause * nother
      allocate(values(nrest), covariance(nrest, nrest), response_idx(nrest), reg_idx(nrest))
      idx = 0
      do eq = 1, model%k
         if (cause(eq)) cycle
         do lag = 1, model%p
            do cause_var = 1, model%k
               if (.not. cause(cause_var)) cycle
               idx = idx + 1
               response_idx(idx) = eq
               reg_idx(idx) = (lag - 1) * model%k + cause_var
               values(idx) = model%coef(eq, reg_idx(idx))
            end do
         end do
      end do
      do b = 1, nrest
         do a = 1, nrest
            covariance(a, b) = sigma(response_idx(a), response_idx(b)) * &
               xtxinv(reg_idx(a), reg_idx(b))
         end do
      end do
      call inverse_matrix(covariance, invcov, invinfo)
      if (invinfo /= 0) then
         info = vars_singular
         return
      end if
      result%statistic = dot_product(values, matmul(invcov, values)) / real(nrest, dp)
      result%df1 = real(nrest, dp)
      result%df2 = real(model%k * model%nobs - model%k * model%nreg, dp)
      result%p_value = r_pf(result%statistic, result%df1, result%df2, lower_tail = .false.)
      info = vars_success
   end subroutine granger_causality

   subroutine instantaneous_causality(model, cause, result, info)
      type(var_model), intent(in) :: model
      logical, intent(in) :: cause(:)
      type(vars_test_result), intent(out) :: result
      integer, intent(out) :: info
      real(dp), allocatable :: sigma(:, :), sigvech(:), cmat(:, :), dmat(:, :), dinv(:, :)
      real(dp), allocatable :: dtdinv(:, :), kron_sigma(:, :), middle(:, :), invmiddle(:, :)
      real(dp), allocatable :: dt(:, :), selector_values(:)
      integer :: k, q, nrest, i, j, row, idx, invinfo
      real(dp) :: denom

      k = model%k
      if (size(cause) /= k .or. count(cause) < 1 .or. count(cause) >= k) then
         info = vars_invalid_argument
         return
      end if
      denom = real(model%nobs - model%nreg, dp)
      allocate(sigma(k, k))
      sigma = matmul(transpose(model%resid), model%resid) / denom
      q = k * (k + 1) / 2
      allocate(sigvech(q))
      call vech_lower(sigma, sigvech)
      nrest = count(cause) * (k - count(cause))
      allocate(cmat(nrest, q), selector_values(nrest))
      cmat = 0.0_dp
      row = 0
      idx = 0
      do j = 1, k
         do i = j, k
            idx = idx + 1
            if (cause(i) .neqv. cause(j)) then
               row = row + 1
               cmat(row, idx) = 1.0_dp
               selector_values(row) = sigma(i, j)
            end if
         end do
      end do
      call duplication_matrix(k, dmat)
      dt = transpose(dmat)
      call inverse_matrix(matmul(dt, dmat), dtdinv, invinfo)
      if (invinfo /= 0) then
         info = vars_singular
         return
      end if
      dinv = matmul(dtdinv, dt)
      kron_sigma = kronecker_product(sigma, sigma)
      middle = 2.0_dp * matmul(cmat, matmul(dinv, matmul(kron_sigma, matmul(transpose(dinv), transpose(cmat)))))
      call inverse_matrix(middle, invmiddle, invinfo)
      if (invinfo /= 0) then
         info = vars_singular
         return
      end if
      result%statistic = real(model%nobs, dp) * dot_product(selector_values, matmul(invmiddle, selector_values))
      result%df1 = real(nrest, dp)
      result%df2 = 0.0_dp
      result%p_value = r_pchisq(result%statistic, result%df1, lower_tail = .false.)
      info = vars_success
   end subroutine instantaneous_causality

end module vars_causality
