! SPDX-License-Identifier: MIT
! Derived from etrm 1.0.2, Copyright (c) 2021 etrm authors.
module etrm_msfc
   use etrm_kinds, only : dp
   use etrm_status, only : etrm_ok, etrm_err_size, etrm_err_argument, &
      etrm_err_linear_solve, set_status
   use etrm_types, only : msfc_result
   implicit none
   private

   public :: msfc, maximum_smoothness_forward_curve

   interface msfc
      module procedure msfc_no_prior
      module procedure msfc_scalar_prior
      module procedure msfc_vector_prior
   end interface msfc

   interface
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import dp
         integer, intent(in) :: n, nrhs, lda, ldb
         integer, intent(out) :: ipiv(*)
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         integer, intent(out) :: info
      end subroutine dgesv
   end interface

contains

   subroutine msfc_no_prior(include, contract, start_day, end_day, price, result, status, message)
      logical, intent(in) :: include(:)
      character(len=*), intent(in) :: contract(:)
      integer, intent(in) :: start_day(:), end_day(:)
      real(dp), intent(in) :: price(:)
      type(msfc_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp) :: prior(1)

      prior(1) = 0.0_dp
      call maximum_smoothness_forward_curve(include, contract, start_day, end_day, price, &
         prior, result, status, message)
   end subroutine msfc_no_prior

   subroutine msfc_scalar_prior(include, contract, start_day, end_day, price, prior, &
      result, status, message)
      logical, intent(in) :: include(:)
      character(len=*), intent(in) :: contract(:)
      integer, intent(in) :: start_day(:), end_day(:)
      real(dp), intent(in) :: price(:), prior
      type(msfc_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp) :: prior_array(1)

      prior_array(1) = prior
      call maximum_smoothness_forward_curve(include, contract, start_day, end_day, price, &
         prior_array, result, status, message)
   end subroutine msfc_scalar_prior

   subroutine msfc_vector_prior(include, contract, start_day, end_day, price, prior, &
      result, status, message)
      logical, intent(in) :: include(:)
      character(len=*), intent(in) :: contract(:)
      integer, intent(in) :: start_day(:), end_day(:)
      real(dp), intent(in) :: price(:), prior(:)
      type(msfc_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message

      call maximum_smoothness_forward_curve(include, contract, start_day, end_day, price, &
         prior, result, status, message)
   end subroutine msfc_vector_prior

   subroutine maximum_smoothness_forward_curve(include, contract, start_day, end_day, &
      price, prior_input, result, status, message)
      logical, intent(in) :: include(:)
      character(len=*), intent(in) :: contract(:)
      integer, intent(in) :: start_day(:), end_day(:)
      real(dp), intent(in) :: price(:), prior_input(:)
      type(msfc_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), allocatable :: h(:, :), a(:, :), b(:), kkt(:, :), rhs(:, :)
      real(dp), allocatable :: tcs(:), tce(:), duration(:)
      integer, allocatable :: knot_days(:), ipiv(:)
      integer :: n_all, m, n, nk, qrows, pcols, total, info
      integer :: i, j, row, left_col, max_end, max_name, seg
      real(dp) :: delta, xk, t1, t2, prior_mean, tv, tvo, dtv

      if (.not. validate_msfc_inputs(include, contract, start_day, end_day, price, &
         prior_input, status, message)) return

      n_all = size(price)
      m = count(include)
      max_end = maxval(end_day, mask=include)
      result%n_days = max_end + 1
      result%n_contracts = m
      allocate(result%day(result%n_days), result%prior(result%n_days), result%curve(result%n_days))
      result%day = [(i - 1, i=1, result%n_days)]
      if (size(prior_input) == 1) then
         result%prior = prior_input(1)
      else
         result%prior = prior_input(1:result%n_days)
      end if

      max_name = 1
      do i = 1, n_all
         if (include(i)) max_name = max(max_name, len_trim(contract(i)))
      end do
      allocate(result%original_index(m), result%start_day(m), result%end_day(m), &
         result%market_price(m), result%computed_price(m))
      allocate(character(len=max_name) :: result%contract(m))
      j = 0
      do i = 1, n_all
         if (.not. include(i)) cycle
         j = j + 1
         result%original_index(j) = i
         result%start_day(j) = start_day(i)
         result%end_day(j) = end_day(i)
         result%market_price(j) = price(i)
         result%contract(j) = trim(contract(i))
      end do

      call make_knots(result%start_day, result%end_day, knot_days)
      nk = size(knot_days)
      n = nk - 1
      result%n_polynomials = n
      allocate(result%knots(nk), result%coefficients(5, n))
      result%knots = real(knot_days, dp) / 365.0_dp

      allocate(tcs(m), tce(m), duration(m))
      tcs = real(result%start_day, dp) / 365.0_dp
      tce = real(result%end_day, dp) / 365.0_dp
      duration = tce - tcs

      pcols = 5 * n
      qrows = 3 * n + m - 2
      total = pcols + qrows
      allocate(h(pcols, pcols), a(qrows, pcols), b(qrows))
      h = 0.0_dp
      a = 0.0_dp
      b = 0.0_dp

      do j = 1, n
         delta = result%knots(j + 1) - result%knots(j)
         left_col = 5 * (j - 1) + 1
         h(left_col, left_col) = 28.8_dp * delta**5
         h(left_col, left_col + 1) = 18.0_dp * delta**4
         h(left_col, left_col + 2) = 8.0_dp * delta**3
         h(left_col + 1, left_col) = h(left_col, left_col + 1)
         h(left_col + 1, left_col + 1) = 12.0_dp * delta**3
         h(left_col + 1, left_col + 2) = 6.0_dp * delta**2
         h(left_col + 2, left_col) = h(left_col, left_col + 2)
         h(left_col + 2, left_col + 1) = h(left_col + 1, left_col + 2)
         h(left_col + 2, left_col + 2) = 4.0_dp * delta
      end do

      row = 1
      do j = 2, n
         xk = result%knots(j)
         left_col = 5 * (j - 2) + 1
         a(row, left_col:left_col + 9) = [ &
            -xk**4, -xk**3, -xk**2, -xk, -1.0_dp, &
             xk**4,  xk**3,  xk**2,  xk,  1.0_dp]
         a(row + 1, left_col:left_col + 9) = [ &
            -4.0_dp*xk**3, -3.0_dp*xk**2, -2.0_dp*xk, -1.0_dp, 0.0_dp, &
             4.0_dp*xk**3,  3.0_dp*xk**2,  2.0_dp*xk,  1.0_dp, 0.0_dp]
         a(row + 2, left_col:left_col + 9) = [ &
            -12.0_dp*xk**2, -6.0_dp*xk, -2.0_dp, 0.0_dp, 0.0_dp, &
             12.0_dp*xk**2,  6.0_dp*xk,  2.0_dp, 0.0_dp, 0.0_dp]
         row = row + 3
      end do

      row = 3 * (n - 1) + 1
      xk = result%knots(n + 1)
      left_col = 5 * (n - 1) + 1
      a(row, left_col:left_col + 4) = [4.0_dp*xk**3, 3.0_dp*xk**2, &
         2.0_dp*xk, 1.0_dp, 0.0_dp]

      row = row + 1
      do i = 1, m
         do j = 1, n
            if (tce(i) > result%knots(j) .and. tcs(i) < result%knots(j + 1)) then
               t1 = max(tcs(i), result%knots(j))
               t2 = min(tce(i), result%knots(j + 1))
               left_col = 5 * (j - 1) + 1
               a(row, left_col:left_col + 4) = [ &
                  (t2**5 - t1**5) / 5.0_dp, (t2**4 - t1**4) / 4.0_dp, &
                  (t2**3 - t1**3) / 3.0_dp, (t2**2 - t1**2) / 2.0_dp, t2 - t1]
            end if
         end do
         prior_mean = sum(result%prior(result%start_day(i) + 1:result%end_day(i) + 1)) / &
            real(result%end_day(i) - result%start_day(i) + 1, dp)
         b(row) = (result%market_price(i) - prior_mean) * duration(i)
         row = row + 1
      end do

      allocate(kkt(total, total), rhs(total, 1), ipiv(total))
      kkt = 0.0_dp
      rhs = 0.0_dp
      kkt(1:pcols, 1:pcols) = 2.0_dp * h
      kkt(1:pcols, pcols + 1:total) = transpose(a)
      kkt(pcols + 1:total, 1:pcols) = a
      rhs(pcols + 1:total, 1) = b
      call dgesv(total, 1, kkt, total, ipiv, rhs, total, info)
      if (info /= 0) then
         call set_status(status, message, etrm_err_linear_solve, &
            "MSFC KKT system is singular or could not be solved")
         return
      end if

      do j = 1, n
         result%coefficients(:, j) = rhs(5*(j - 1) + 1:5*j, 1)
      end do

      dtv = 0.0001_dp / 365.0_dp
      do i = 1, result%n_days
         tv = real(result%day(i), dp) / 365.0_dp
         tvo = tv + dtv
         seg = segment_for_time(tv, result%knots)
         result%curve(i) = polynomial_integral(result%coefficients(:, seg), tv, tvo) / dtv + &
            result%prior(i)
      end do

      do i = 1, m
         result%computed_price(i) = 0.0_dp
         do j = 1, n
            t1 = max(tcs(i), result%knots(j))
            t2 = min(tce(i), result%knots(j + 1))
            if (t2 > t1) then
               result%computed_price(i) = result%computed_price(i) + &
                  polynomial_integral(result%coefficients(:, j), t1, t2)
            end if
         end do
         prior_mean = sum(result%prior(result%start_day(i) + 1:result%end_day(i) + 1)) / &
            real(result%end_day(i) - result%start_day(i) + 1, dp)
         result%computed_price(i) = result%computed_price(i) / duration(i) + prior_mean
      end do

      call set_status(status, message, etrm_ok, "")
   end subroutine maximum_smoothness_forward_curve

   logical function validate_msfc_inputs(include, contract, start_day, end_day, price, &
      prior_input, status, message) result(valid)
      logical, intent(in) :: include(:)
      character(len=*), intent(in) :: contract(:)
      integer, intent(in) :: start_day(:), end_day(:)
      real(dp), intent(in) :: price(:), prior_input(:)
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      integer :: n, max_end

      valid = .false.
      n = size(price)
      if (n < 1) then
         call set_status(status, message, etrm_err_size, "At least one contract is required")
      else if (size(include) /= n .or. size(contract) /= n .or. size(start_day) /= n .or. &
         size(end_day) /= n) then
         call set_status(status, message, etrm_err_size, &
            "include, contract, start_day, end_day, and price must have equal lengths")
      else if (count(include) < 1) then
         call set_status(status, message, etrm_err_argument, "At least one contract must be included")
      else if (any((start_day < 0) .and. include)) then
         call set_status(status, message, etrm_err_argument, &
            "Included contracts cannot start before the curve trade date")
      else if (any(((end_day - start_day) < 1) .and. include)) then
         call set_status(status, message, etrm_err_argument, &
            "Each included contract must have end_day at least one day after start_day")
      else
         max_end = maxval(end_day, mask=include)
         if (size(prior_input) /= 1 .and. size(prior_input) < max_end + 1) then
            call set_status(status, message, etrm_err_size, &
               "A nonconstant prior must cover day zero through the final included end day")
         else
            valid = .true.
         end if
      end if
   end function validate_msfc_inputs

   subroutine make_knots(start_day, end_day, knots)
      integer, intent(in) :: start_day(:), end_day(:)
      integer, allocatable, intent(out) :: knots(:)
      integer, allocatable :: work(:)
      integer :: i, j, key, n_unique

      allocate(work(1 + 2 * size(start_day)))
      work(1) = 0
      do i = 1, size(start_day)
         work(2*i) = start_day(i)
         work(2*i + 1) = end_day(i)
      end do
      do i = 2, size(work)
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      n_unique = 1
      do i = 2, size(work)
         if (work(i) /= work(n_unique)) then
            n_unique = n_unique + 1
            work(n_unique) = work(i)
         end if
      end do
      allocate(knots(n_unique))
      knots = work(1:n_unique)
   end subroutine make_knots

   pure integer function segment_for_time(t, knots) result(segment)
      real(dp), intent(in) :: t, knots(:)
      integer :: j

      segment = size(knots) - 1
      do j = 1, size(knots) - 1
         if (t <= knots(j + 1) + 32.0_dp * epsilon(t)) then
            segment = j
            exit
         end if
      end do
   end function segment_for_time

   pure real(dp) function polynomial_integral(coef, x1, x2) result(value)
      real(dp), intent(in) :: coef(5), x1, x2
      value = coef(1) / 5.0_dp * (x2**5 - x1**5) + &
         coef(2) / 4.0_dp * (x2**4 - x1**4) + &
         coef(3) / 3.0_dp * (x2**3 - x1**3) + &
         coef(4) / 2.0_dp * (x2**2 - x1**2) + &
         coef(5) * (x2 - x1)
   end function polynomial_integral

end module etrm_msfc
