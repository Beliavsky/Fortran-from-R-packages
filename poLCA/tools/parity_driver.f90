program parity_driver
   use polca, only : dp, polca_model, polca_fit, polca_update_prior, polca_rmulti
   implicit none

   call nocov_case()
   call regression_case()

contains

   subroutine nocov_case()
      integer, parameter :: n = 2400
      integer :: y(n, 3), k(3), cell(3), i, idx, rep
      real(dp) :: start(2, 2, 3), prob, share(2)
      type(polca_model) :: model

      k = 2
      start(1, :, 1) = [0.90_dp, 0.10_dp]
      start(1, :, 2) = [0.80_dp, 0.20_dp]
      start(1, :, 3) = [0.85_dp, 0.15_dp]
      start(2, :, 1) = [0.15_dp, 0.85_dp]
      start(2, :, 2) = [0.25_dp, 0.75_dp]
      start(2, :, 3) = [0.20_dp, 0.80_dp]
      share = [0.4_dp, 0.6_dp]
      idx = 0
      do i = 0, 7
         cell = [1 + iand(i, 1), 1 + iand(ishft(i, -1), 1), 1 + iand(ishft(i, -2), 1)]
         prob = mixture_cell(start, share, cell)
         do rep = 1, nint(prob * real(n, dp))
            if (idx < n) then
               idx = idx + 1
               y(idx, :) = cell
            end if
         end do
      end do
      do while (idx < n)
         idx = idx + 1
         y(idx, :) = [2, 2, 2]
      end do
      call polca_fit(y, k, 2, model, probs_start=start, maxiter=2000, tol=1.0e-13_dp, calc_se=.true.)
      write(*, '(a,1x,es24.16)') 'nocov_loglik', model%loglik
      write(*, '(a,2(1x,es24.16))') 'nocov_share', model%class_share
      write(*, '(a,3(1x,es24.16))') 'nocov_p2_c1', model%probs(1, 2, :)
      write(*, '(a,3(1x,es24.16))') 'nocov_p2_c2', model%probs(2, 2, :)
      write(*, '(a,2(1x,es24.16))') 'nocov_share_se', model%class_share_se
      write(*, '(a,3(1x,es24.16))') 'nocov_p2_se_c1', model%probs_se(1, 2, :)
      write(*, '(a,3(1x,es24.16))') 'nocov_p2_se_c2', model%probs_se(2, 2, :)
   end subroutine nocov_case

   subroutine regression_case()
      integer, parameter :: n = 1200
      integer :: y(n, 3), k(3), i
      real(dp) :: probs(2, 2, 3), x(n, 2), beta(2, 1), prior(n, 2)
      real(dp) :: u(n), pmat(n, 2), start(2, 2, 3)
      integer :: cls(n)
      type(polca_model) :: model

      k = 2
      probs(1, :, 1) = [0.90_dp, 0.10_dp]
      probs(1, :, 2) = [0.85_dp, 0.15_dp]
      probs(1, :, 3) = [0.80_dp, 0.20_dp]
      probs(2, :, 1) = [0.10_dp, 0.90_dp]
      probs(2, :, 2) = [0.20_dp, 0.80_dp]
      probs(2, :, 3) = [0.15_dp, 0.85_dp]
      beta(:, 1) = [-0.25_dp, 1.1_dp]
      do i = 1, n
         x(i, 1) = 1.0_dp
         x(i, 2) = -1.5_dp + 3.0_dp * real(i - 1, dp) / real(n - 1, dp)
      end do
      call polca_update_prior(beta, x, prior)
      call deterministic_uniform(n, 17, u)
      call polca_rmulti(prior, u, cls)
      do i = 1, 3
         pmat(:, 1) = probs(cls, 1, i)
         pmat(:, 2) = probs(cls, 2, i)
         call deterministic_uniform(n, 40 + i, u)
         call polca_rmulti(pmat, u, y(:, i))
      end do
      start = probs
      call polca_fit(y, k, 2, model, x=x, probs_start=start, maxiter=2000, tol=1.0e-12_dp, calc_se=.true.)
      write(*, '(a,1x,es24.16)') 'reg_loglik', model%loglik
      write(*, '(a,2(1x,es24.16))') 'reg_beta', model%coeff(:, 1)
      write(*, '(a,3(1x,es24.16))') 'reg_p2_c1', model%probs(1, 2, :)
      write(*, '(a,3(1x,es24.16))') 'reg_p2_c2', model%probs(2, 2, :)
      write(*, '(a,2(1x,es24.16))') 'reg_beta_se', model%coeff_se(:, 1)
      write(*, '(a,2(1x,es24.16))') 'reg_share_se', model%class_share_se
   end subroutine regression_case

   pure real(dp) function mixture_cell(probs, share, cell) result(p)
      real(dp), intent(in) :: probs(:, :, :) !! Class-by-category-by-item probabilities for the parity mixture.
      real(dp), intent(in) :: share(:) !! Class probabilities for the parity mixture.
      integer, intent(in) :: cell(:) !! Complete response pattern to evaluate.
      real(dp) :: pr
      integer :: j, r

      p = 0.0_dp
      do r = 1, size(share)
         pr = share(r)
         do j = 1, size(cell)
            pr = pr * probs(r, cell(j), j)
         end do
         p = p + pr
      end do
   end function mixture_cell

   pure subroutine deterministic_uniform(n, offset, u)
      use, intrinsic :: iso_fortran_env, only : int64
      integer, intent(in) :: n !! Number of deterministic pseudo-uniform values.
      integer, intent(in) :: offset !! Stream selector for the deterministic sequence.
      real(dp), intent(out) :: u(:) !! Deterministic values in the open unit interval.
      integer(int64) :: state
      integer :: i

      state = modulo(123457_int64 + 104729_int64 * int(offset, int64), 2147483646_int64) + 1_int64
      do i = 1, n
         state = modulo(48271_int64 * state, 2147483647_int64)
         u(i) = real(state, dp) / 2147483647.0_dp
      end do
   end subroutine deterministic_uniform

end program parity_driver
