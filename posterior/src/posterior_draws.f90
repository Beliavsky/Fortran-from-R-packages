! SPDX-License-Identifier: BSD-3-Clause
! Numeric translation of posterior 1.7.0. See NOTICE and upstream provenance.
module posterior_draws
   use r_kinds, only : dp
   use r_descriptive, only : r_mean, r_variance, r_sd
   use r_stability, only : r_log_sum_exp
   use r_ordering, only : r_order, r_sort_values_in_place
   use r_quantiles, only : r_quantile_type7
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   implicit none
   private
   public :: dp, draws_to_matrix, matrix_to_draws, merge_draws_chains, split_draws_chains
   public :: bind_draws_arrays, subset_draws_array, thin_draws_array, order_draws_array
   public :: resample_draws_array, extract_draws_variable, extract_draws_variable_matrix
   public :: rvar_summary, rvar_quantiles, rvar_ranges, rvar_expectation, rvar_probability
   public :: weight_draws_array, extract_draws_weights, rvar_over_draws_summary
   public :: rvar_cdf, rvar_matmul, rvar_ifelse, rvar_special_masks, rvar_constant
contains

   !> Flatten iterations and chains, preserving R's chain-major draw order.
   pure subroutine draws_to_matrix(x, values)
      real(dp), intent(in) :: x(:, :, :) !! Draws with dimensions (iteration, chain, variable).
      real(dp), allocatable, intent(out) :: values(:, :) !! Flattened (draw, variable) values.
      values = reshape(x, [size(x, 1) * size(x, 2), size(x, 3)])
   end subroutine draws_to_matrix

   !> Restore the chain dimension of flattened numeric draws.
   pure subroutine matrix_to_draws(x, nchains, values, status)
      real(dp), intent(in) :: x(:, :) !! Flattened (draw, variable) values, chains in contiguous blocks.
      integer, intent(in) :: nchains !! Positive number of chains, dividing the number of draws.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! (iteration, chain, variable) values, empty on error.
      integer, intent(out), optional :: status !! Zero on success, one for incompatible chain count.
      if (present(status)) status = 1
      allocate(values(0, 0, 0))
      if (nchains < 1) return
      if (mod(size(x, 1), nchains) /= 0) return
      values = reshape(x, [size(x, 1) / nchains, nchains, size(x, 2)])
      if (present(status)) status = 0
   end subroutine matrix_to_draws

   pure subroutine merge_draws_chains(x, values)
      real(dp), intent(in) :: x(:, :, :) !! Draws with dimensions (iteration, chain, variable).
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Draws with one chain and original draw order.
      values = reshape(x, [size(x, 1) * size(x, 2), 1, size(x, 3)])
   end subroutine merge_draws_chains

   !> Public split_chains semantics: discard the last iteration when odd.
   pure subroutine split_draws_chains(x, values)
      real(dp), intent(in) :: x(:, :, :) !! Draws with dimensions (iteration, chain, variable).
      real(dp), allocatable, intent(out) :: values(:, :, :) !! First halves followed by second halves as chains.
      integer :: half, chains
      half = size(x, 1) / 2
      chains = size(x, 2)
      allocate(values(half, 2 * chains, size(x, 3)))
      values(:, :chains, :) = x(:half, :, :)
      values(:, chains + 1:, :) = x(half + 1:2 * half, :, :)
   end subroutine split_draws_chains

   pure subroutine bind_draws_arrays(x, y, along, values, status)
      real(dp), intent(in) :: x(:, :, :) !! First (iteration, chain, variable) draws array.
      real(dp), intent(in) :: y(:, :, :) !! Second draws array, equal in the two non-binding dimensions.
      character(len=*), intent(in) :: along !! Binding dimension: iteration, chain, or variable.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Concatenated values, or a zero-sized array on error.
      integer, intent(out), optional :: status !! Zero on success, one for an unknown dimension or incompatible sizes.
      integer :: nx(3), ny(3)
      nx = shape(x)
      ny = shape(y)
      allocate(values(0, 0, 0))
      if (present(status)) status = 1
      select case (along)
      case ('iteration')
         if (any(nx(2:3) /= ny(2:3))) return
         deallocate(values)
         allocate(values(nx(1) + ny(1), nx(2), nx(3)))
         values(:nx(1), :, :) = x
         values(nx(1) + 1:, :, :) = y
      case ('chain')
         if (nx(1) /= ny(1) .or. nx(3) /= ny(3)) return
         deallocate(values)
         allocate(values(nx(1), nx(2) + ny(2), nx(3)))
         values(:, :nx(2), :) = x
         values(:, nx(2) + 1:, :) = y
      case ('variable')
         if (any(nx(1:2) /= ny(1:2))) return
         deallocate(values)
         allocate(values(nx(1), nx(2), nx(3) + ny(3)))
         values(:, :, :nx(3)) = x
         values(:, :, nx(3) + 1:) = y
      case default
         return
      end select
      if (present(status)) status = 0
   end subroutine bind_draws_arrays

   !> Positive positional selections, with sorting of iteration/chain/draw indices as in R.
   pure subroutine subset_draws_array(x, values, iteration, chain, variable, draw, unique, exclude, status)
      real(dp), intent(in) :: x(:, :, :) !! Draws in canonical (iteration, chain, variable) order.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Selected values, or a zero-sized array on invalid indices.
      integer, intent(in), optional :: iteration(:) !! One-based iterations, sorted before selection; absent selects all.
      integer, intent(in), optional :: chain(:) !! One-based chains, sorted before selection; absent selects all.
      integer, intent(in), optional :: variable(:) !! Resolved one-based scalar variable positions, preserving supplied order.
      integer, intent(in), optional :: draw(:) !! One-based pooled draws; cannot be combined with chain or iteration.
      logical, intent(in), optional :: unique !! Remove repeated iteration, chain, and draw selections; default true.
      logical, intent(in), optional :: exclude !! Select the complement of each supplied index vector; default false.
      integer, intent(out), optional :: status !! Zero on success, one for out-of-range or conflicting indices.
      integer, allocatable :: ii(:), cc(:), vv(:)
      real(dp), allocatable :: pooled(:, :, :)
      logical :: distinct, complement, ok
      distinct = .true.
      if (present(unique)) distinct = unique
      complement = .false.
      if (present(exclude)) complement = exclude
      if (present(status)) status = 1
      allocate(values(0, 0, 0))
      if (present(draw)) then
         if (present(iteration) .or. present(chain)) return
      end if
      call select_positions(size(x, 3), variable, .true., complement, .false., vv, ok)
      if (.not. ok) return
      if (present(draw)) then
         call select_positions(size(x, 1) * size(x, 2), draw, distinct, complement, .true., ii, ok)
         if (.not. ok) return
         call merge_draws_chains(x, pooled)
         values = pooled(ii, :, vv)
      else
         call select_positions(size(x, 1), iteration, distinct, complement, .true., ii, ok)
         if (.not. ok) return
         call select_positions(size(x, 2), chain, distinct, complement, .true., cc, ok)
         if (.not. ok) return
         values = x(ii, cc, vv)
      end if
      if (present(status)) status = 0
   end subroutine subset_draws_array

   pure subroutine select_positions(n, requested, unique, exclude, sorted, selected, ok)
      integer, intent(in) :: n !! Length of the dimension being selected.
      integer, intent(in), optional :: requested(:) !! Requested one-based positions; absent means all positions.
      logical, intent(in) :: unique !! Whether to discard repeated requested positions.
      logical, intent(in) :: exclude !! Whether to select the complement of supplied positions.
      logical, intent(in) :: sorted !! Whether to sort requested positions into ascending order.
      integer, allocatable, intent(out) :: selected(:) !! Validated selection, empty on error.
      logical, intent(out) :: ok !! True if every requested position is in bounds.
      integer, allocatable :: ids(:), ordering(:)
      logical, allocatable :: seen(:)
      integer :: i, k
      ok = .false.
      allocate(selected(0))
      if (.not. present(requested)) then
         selected = [(i, i = 1, n)]
         ok = .true.
         return
      end if
      if (any(requested < 1) .or. any(requested > n)) return
      allocate(seen(n))
      seen = .false.
      allocate(ids(size(requested)))
      k = 0
      do i = 1, size(requested)
         if (unique .and. seen(requested(i))) cycle
         k = k + 1
         ids(k) = requested(i)
         seen(requested(i)) = .true.
      end do
      if (exclude) then
         selected = pack([(i, i = 1, n)], .not. seen)
      else if (sorted) then
         call r_order(real(ids(:k), dp), ordering)
         selected = ids(ordering)
      else
         selected = ids(:k)
      end if
      ok = .true.
   end subroutine select_positions

   pure subroutine thin_draws_array(x, thin, values, status)
      real(dp), intent(in) :: x(:, :, :) !! Draws with dimensions (iteration, chain, variable).
      real(dp), intent(in) :: thin !! Finite spacing at least one and at most the number of iterations, unless exactly one.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Draws at round(seq(1,niterations,by=thin)), with ties to even.
      integer, intent(out), optional :: status !! Zero on success, one for invalid spacing; no automatic ESS-based spacing.
      integer, allocatable :: ids(:)
      integer :: i, n
      allocate(values(0, 0, 0))
      if (present(status)) status = 1
      if (.not. ieee_is_finite(thin)) return
      if (thin == 1.0_dp) then
         values = x
      else
         if (thin < 1.0_dp .or. thin > real(size(x, 1), dp)) return
         n = int(floor(real(size(x, 1) - 1, dp) / thin + 1.0e-10_dp)) + 1
         allocate(ids(n))
         do i = 1, n
            ids(i) = round_even(min(real(size(x, 1), dp), 1.0_dp + real(i - 1, dp) * thin))
         end do
         values = x(ids, :, :)
      end if
      if (present(status)) status = 0
   end subroutine thin_draws_array

   pure integer function round_even(x) result(value)
      real(dp), intent(in) :: x !! Nonnegative finite number in default integer range.
      real(dp) :: fraction
      value = int(floor(x))
      fraction = x - real(value, dp)
      if (fraction > 0.5_dp) value = value + 1
      if (fraction == 0.5_dp .and. mod(value, 2) /= 0) value = value + 1
   end function round_even

   pure subroutine order_draws_array(x, iteration_ids, chain_ids, values, status)
      real(dp), intent(in) :: x(:, :, :) !! Possibly out-of-order (iteration, chain, variable) draws.
      integer, intent(in) :: iteration_ids(:) !! Iteration labels, one for each first-dimension position.
      integer, intent(in) :: chain_ids(:) !! Chain labels, one for each second-dimension position.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Draws ordered stably by ascending iteration and chain labels.
      integer, intent(out), optional :: status !! Zero on success, one if label counts do not match dimensions.
      integer, allocatable :: ii(:), cc(:)
      allocate(values(0, 0, 0))
      if (present(status)) status = 1
      if (size(iteration_ids) /= size(x, 1) .or. size(chain_ids) /= size(x, 2)) return
      call r_order(real(iteration_ids, dp), ii)
      call r_order(real(chain_ids, dp), cc)
      values = x(ii, cc, :)
      if (present(status)) status = 0
   end subroutine order_draws_array

   !> Numeric resampling using posterior's deterministic or stratified algorithm.
   pure subroutine resample_draws_array(x, weights, ndraws, values, method, uniforms, status)
      real(dp), intent(in) :: x(:, :, :) !! Draws with dimensions (iteration, chain, variable).
      real(dp), intent(in) :: weights(:) !! Finite nonnegative weights in pooled draw order, with positive sum.
      integer, intent(in) :: ndraws !! Nonnegative number of output draws.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Resampled (ndraws, 1, variable) array, sorted as in subset_draws.
      character(len=*), intent(in), optional :: method !! 'deterministic' (default) or 'stratified'.
      real(dp), intent(in), optional :: uniforms(:) !! Stratified U[0,1] sequence in consumption order, at least size(weights).
      integer, intent(out), optional :: status !! Zero on success, one for invalid arguments or numerical failure.
      real(dp), allocatable :: w(:)
      integer, allocatable :: ids(:)
      real(dp) :: carry, total, fraction
      integer :: i, j, a, iu, inner_status
      character(len=:), allocatable :: selected_method
      allocate(values(0, 0, 0))
      if (present(status)) status = 1
      if (ndraws < 0 .or. size(weights) /= size(x, 1) * size(x, 2)) return
      if (size(weights) == 0) return
      if (any(.not. ieee_is_finite(weights)) .or. any(weights < 0.0_dp)) return
      total = sum(weights)
      if (.not. ieee_is_finite(total) .or. total <= 0.0_dp) return
      selected_method = 'deterministic'
      if (present(method)) selected_method = method
      if (selected_method /= 'deterministic' .and. selected_method /= 'stratified') return
      if (selected_method == 'stratified') then
         if (.not. present(uniforms)) return
         if (size(uniforms) < size(weights)) return
         if (any(.not. ieee_is_finite(uniforms))) return
         if (any(uniforms < 0.0_dp) .or. any(uniforms > 1.0_dp)) return
      end if
      w = real(ndraws, dp) * (weights / total)
      allocate(ids(ndraws))
      j = 0
      iu = 0
      if (selected_method == 'deterministic') then
         carry = 0.5_dp
         do i = 1, size(w)
            a = int(floor(w(i)))
            if (j + a > ndraws) return
            ids(j + 1:j + a) = i
            j = j + a
            carry = carry + w(i) - real(a, dp)
            if (carry >= 1.0_dp) then
               if (j == ndraws) return
               j = j + 1
               ids(j) = i
               carry = carry - 1.0_dp
            end if
         end do
      else
         carry = 0.0_dp
         do i = 1, size(w)
            carry = carry + w(i)
            if (carry >= 1.0_dp) then
               a = int(floor(carry))
               if (j + a > ndraws) return
               carry = carry - real(a, dp)
               ids(j + 1:j + a) = i
               j = j + a
            end if
            if (j < ndraws) then
               iu = iu + 1
               fraction = uniforms(iu)
               if (carry >= fraction) then
                  carry = carry - 1.0_dp
                  j = j + 1
                  ids(j) = i
               end if
            end if
         end do
      end if
      if (j /= ndraws) return
      call subset_draws_array(x, values, draw = ids, unique = .false., status = inner_status)
      if (present(status)) status = inner_status
   end subroutine resample_draws_array

   pure subroutine extract_draws_variable(x, variable, values, status)
      real(dp), intent(in) :: x(:, :, :) !! Draws with dimensions (iteration, chain, variable).
      integer, intent(in) :: variable !! One-based position of the scalar variable to extract.
      real(dp), allocatable, intent(out) :: values(:) !! Pooled draw vector in chain-major order, empty on error.
      integer, intent(out), optional :: status !! Zero on success, one when the variable index is out of bounds.
      allocate(values(0))
      if (present(status)) status = 1
      if (variable < 1 .or. variable > size(x, 3)) return
      values = reshape(x(:, :, variable), [size(x, 1) * size(x, 2)])
      if (present(status)) status = 0
   end subroutine extract_draws_variable

   pure subroutine extract_draws_variable_matrix(x, variable, values, status)
      real(dp), intent(in) :: x(:, :, :) !! Draws with dimensions (iteration, chain, variable).
      integer, intent(in) :: variable !! One-based position of the scalar variable to extract.
      real(dp), allocatable, intent(out) :: values(:, :) !! (iteration, chain) values, empty on error.
      integer, intent(out), optional :: status !! Zero on success, one when the variable index is out of bounds.
      allocate(values(0, 0))
      if (present(status)) status = 1
      if (variable < 1 .or. variable > size(x, 3)) return
      values = x(:, :, variable)
      if (present(status)) status = 0
   end subroutine extract_draws_variable_matrix

   !> Represent numeric constants as an rvar with one draw for broadcasting.
   pure subroutine rvar_constant(x, values)
      real(dp), intent(in) :: x(:) !! Numeric constant elements, in R column-major flattened order.
      real(dp), allocatable, intent(out) :: values(:, :) !! Shape (1, size(x)); one deterministic draw per element.
      values = reshape(x, [1, size(x)])
   end subroutine rvar_constant

   !> Apply a reduction across elements independently within each draw.
   pure subroutine rvar_summary(x, method, values, na_rm, constant, status)
      real(dp), intent(in) :: x(:, :) !! Numeric rvar values with dimensions (draw, element); missing values are IEEE NaNs.
      character(len=*), intent(in) :: method !! mean, median, sum, prod, min, max, sd, var, mad, all, or any.
      real(dp), allocatable, intent(out) :: values(:) !! One summary per draw, with logical results encoded as 0, 1, or NaN.
      logical, intent(in), optional :: na_rm !! Remove missing elements before each reduction; default false.
      real(dp), intent(in), optional :: constant !! MAD multiplier, default 1.4826; ignored by other reductions.
      integer, intent(out), optional :: status !! Zero on success, one for an unsupported method.
      real(dp), allocatable :: row(:)
      real(dp) :: factor, median
      logical :: remove_na, missing
      integer :: i
      allocate(values(size(x, 1)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = 1
      select case (method)
      case ('mean', 'median', 'sum', 'prod', 'min', 'max', 'sd', 'var', 'mad', 'all', 'any')
      case default
         return
      end select
      factor = 1.4826_dp
      if (present(constant)) factor = constant
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      do i = 1, size(x, 1)
         missing = any(ieee_is_nan(x(i, :)))
         if (remove_na) then
            row = pack(x(i, :), .not. ieee_is_nan(x(i, :)))
         else
            row = x(i, :)
         end if
         if (method == 'all') then
            if (any(row == 0.0_dp)) then
               values(i) = 0.0_dp
            else if (.not. missing .or. remove_na) then
               values(i) = 1.0_dp
            end if
            cycle
         else if (method == 'any') then
            if (any(.not. ieee_is_nan(row) .and. row /= 0.0_dp)) then
               values(i) = 1.0_dp
            else if (.not. missing .or. remove_na) then
               values(i) = 0.0_dp
            end if
            cycle
         end if
         if (missing .and. .not. remove_na) cycle
         select case (method)
         case ('mean')
            values(i) = r_mean(row)
         case ('median')
            values(i) = numeric_quantile(row, 0.5_dp)
         case ('sum')
            values(i) = sum(row)
         case ('prod')
            values(i) = product(row)
         case ('min')
            values(i) = ieee_value(0.0_dp, ieee_positive_inf)
            if (size(row) > 0) values(i) = minval(row)
         case ('max')
            values(i) = ieee_value(0.0_dp, ieee_negative_inf)
            if (size(row) > 0) values(i) = maxval(row)
         case ('var')
            values(i) = r_variance(row)
         case ('sd')
            values(i) = r_sd(row)
         case ('mad')
            median = numeric_quantile(row, 0.5_dp)
            values(i) = factor * numeric_quantile(abs(row - median), 0.5_dp)
         end select
      end do
      if (present(status)) status = 0
   end subroutine rvar_summary

   !> Type-7 sample quantiles. Shared core covers finite samples; IEEE endpoints need explicit handling.
   pure real(dp) function numeric_quantile(x, probability) result(value)
      real(dp), intent(in) :: x(:) !! NaN-free samples, including infinities; an empty sample returns NaN.
      real(dp), intent(in) :: probability !! Quantile probability in [0,1], validated by the caller.
      real(dp), allocatable :: ordered(:)
      real(dp) :: h, fraction
      integer :: j, n
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      n = size(x)
      if (n == 0 .or. any(ieee_is_nan(x))) return
      if (all(ieee_is_finite(x))) then
         value = r_quantile_type7(x, probability)
         return
      end if
      ordered = x
      call r_sort_values_in_place(ordered)
      h = 1.0_dp + real(n - 1, dp) * probability
      j = int(floor(h))
      fraction = h - real(j, dp)
      value = ordered(j)
      if (j == n .or. fraction == 0.0_dp) return
      if (ordered(j) == ordered(j + 1)) return
      value = (1.0_dp - fraction) * ordered(j) + fraction * ordered(j + 1)
   end function numeric_quantile

   pure subroutine rvar_quantiles(x, probs, values, na_rm, status)
      real(dp), intent(in) :: x(:, :) !! Numeric rvar samples with dimensions (draw, element).
      real(dp), intent(in) :: probs(:) !! Finite probabilities in [0,1]; numeric type-7 quantiles only.
      real(dp), allocatable, intent(out) :: values(:, :) !! Shape (draw, probability), NaN on invalid input.
      logical, intent(in), optional :: na_rm !! Remove NaNs; default false propagates missingness per draw.
      integer, intent(out), optional :: status !! Zero on success, one for invalid probabilities.
      real(dp), allocatable :: row(:)
      logical :: remove_na
      integer :: i, j
      allocate(values(size(x, 1), size(probs)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = 1
      if (any(.not. ieee_is_finite(probs))) return
      if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp)) return
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      do i = 1, size(x, 1)
         if (.not. remove_na .and. any(ieee_is_nan(x(i, :)))) cycle
         row = pack(x(i, :), .not. ieee_is_nan(x(i, :)))
         do j = 1, size(probs)
            values(i, j) = numeric_quantile(row, probs(j))
         end do
      end do
      if (present(status)) status = 0
   end subroutine rvar_quantiles

   pure subroutine rvar_ranges(x, values, na_rm)
      real(dp), intent(in) :: x(:, :) !! Numeric rvar samples with dimensions (draw, element).
      real(dp), allocatable, intent(out) :: values(:, :) !! Shape (draw, 2), holding minimum and maximum for each draw.
      logical, intent(in), optional :: na_rm !! Remove NaNs before extrema; default false.
      real(dp), allocatable :: lower(:), upper(:)
      call rvar_summary(x, 'min', lower, na_rm)
      call rvar_summary(x, 'max', upper, na_rm)
      allocate(values(size(x, 1), 2))
      values(:, 1) = lower
      values(:, 2) = upper
   end subroutine rvar_ranges

   pure subroutine rvar_expectation(x, values, na_rm)
      real(dp), intent(in) :: x(:, :) !! Numeric rvar samples with dimensions (draw, element).
      real(dp), allocatable, intent(out) :: values(:) !! Mean across draws for each element, NaN for no available draws.
      logical, intent(in), optional :: na_rm !! Remove NaNs independently per element; default false.
      integer :: j
      allocate(values(size(x, 2)))
      do j = 1, size(x, 2)
         values(j) = r_mean(x(:, j), na_rm)
      end do
   end subroutine rvar_expectation

   pure subroutine rvar_probability(x, values, na_rm, status)
      real(dp), intent(in) :: x(:, :) !! Logical draws encoded as 0=false, 1=true, NaN=missing; shape (draw, element).
      real(dp), allocatable, intent(out) :: values(:) !! Empirical probability of true for each element, NaN on error.
      logical, intent(in), optional :: na_rm !! Remove missing logical draws; default false.
      integer, intent(out), optional :: status !! Zero on success, one if a nonmissing input is not zero or one.
      allocate(values(size(x, 2)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = 1
      if (any(.not. ieee_is_nan(x) .and. x /= 0.0_dp .and. x /= 1.0_dp)) return
      call rvar_expectation(x, values, na_rm)
      if (present(status)) status = 0
   end subroutine rvar_probability

   pure subroutine rvar_cdf(x, q, values, status)
      real(dp), intent(in) :: x(:, :) !! Numeric rvar samples (draw, element); missing draws are omitted as in ecdf.
      real(dp), intent(in) :: q(:) !! Evaluation thresholds, including infinities; missing thresholds return NaN.
      real(dp), allocatable, intent(out) :: values(:, :) !! Empirical P(X<=q), shape (threshold, element).
      integer, intent(out), optional :: status !! Zero on success, one if an element has no nonmissing observations.
      real(dp), allocatable :: draws(:)
      integer :: i, j
      allocate(values(size(q), size(x, 2)))
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = 1
      do j = 1, size(x, 2)
         if (count(.not. ieee_is_nan(x(:, j))) == 0) return
      end do
      do j = 1, size(x, 2)
         draws = pack(x(:, j), .not. ieee_is_nan(x(:, j)))
         do i = 1, size(q)
            if (ieee_is_nan(q(i))) cycle
            values(i, j) = real(count(draws <= q(i)), dp) / real(size(draws), dp)
         end do
      end do
      if (present(status)) status = 0
   end subroutine rvar_cdf

   pure subroutine rvar_matmul(x, y, values, status)
      real(dp), intent(in) :: x(:, :, :) !! Left matrix rvar, shape (draw, row, inner); one draw broadcasts as a constant.
      real(dp), intent(in) :: y(:, :, :) !! Right matrix rvar, shape (draw, inner, column); one draw broadcasts as a constant.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Per-draw matrix product, shape (draw, row, column); empty on error.
      integer, intent(out), optional :: status !! Zero on success, one for incompatible matrix or draw dimensions.
      integer :: i, draws, ix, iy
      allocate(values(0, 0, 0))
      if (present(status)) status = 1
      if (size(x, 3) /= size(y, 2)) return
      draws = common_size(size(x, 1), size(y, 1))
      if (draws < 1) return
      deallocate(values)
      allocate(values(draws, size(x, 2), size(y, 3)))
      do i = 1, draws
         ix = min(i, size(x, 1))
         iy = min(i, size(y, 1))
         values(i, :, :) = matmul(x(ix, :, :), y(iy, :, :))
      end do
      if (present(status)) status = 0
   end subroutine rvar_matmul

   pure subroutine rvar_ifelse(test, yes, no, values, status)
      real(dp), intent(in) :: test(:, :) !! Logical draws encoded as 0, 1, or NaN, with dimensions (draw, element).
      real(dp), intent(in) :: yes(:, :) !! Selected values for true draws; each dimension can have length one to broadcast.
      real(dp), intent(in) :: no(:, :) !! Selected values for false draws; each dimension can have length one to broadcast.
      real(dp), allocatable, intent(out) :: values(:, :) !! Selections: missing tests give NaN; invalid shapes give empty output.
      integer, intent(out), optional :: status !! Zero on success, one for incompatible shapes or nonlogical test values.
      integer :: nd, ne, i, j
      real(dp) :: condition
      allocate(values(0, 0))
      if (present(status)) status = 1
      if (any(.not. ieee_is_nan(test) .and. test /= 0.0_dp .and. test /= 1.0_dp)) return
      nd = common_size(common_size(size(test, 1), size(yes, 1)), size(no, 1))
      ne = common_size(common_size(size(test, 2), size(yes, 2)), size(no, 2))
      if (nd < 1 .or. ne < 0) return
      deallocate(values)
      allocate(values(nd, ne))
      do j = 1, ne
         do i = 1, nd
            condition = test(min(i, size(test, 1)), min(j, size(test, 2)))
            if (ieee_is_nan(condition)) then
               values(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
            else if (condition == 1.0_dp) then
               values(i, j) = yes(min(i, size(yes, 1)), min(j, size(yes, 2)))
            else
               values(i, j) = no(min(i, size(no, 1)), min(j, size(no, 2)))
            end if
         end do
      end do
      if (present(status)) status = 0
   end subroutine rvar_ifelse

   pure integer function common_size(nx, ny) result(n)
      integer, intent(in) :: nx !! First dimension length, or -1 to propagate an earlier mismatch.
      integer, intent(in) :: ny !! Second dimension length, or -1 to propagate an earlier mismatch.
      n = -1
      if (nx < 0 .or. ny < 0) return
      if (nx == ny) then
         n = nx
      else if (nx == 1) then
         n = ny
      else if (ny == 1) then
         n = nx
      end if
   end function common_size

   pure subroutine rvar_special_masks(x, finite, infinite, missing)
      real(dp), intent(in) :: x(:, :) !! Numeric rvar draws, including IEEE NaNs and infinities.
      logical, allocatable, intent(out) :: finite(:, :) !! Per-draw finite-value mask with the same shape as x.
      logical, allocatable, intent(out) :: infinite(:, :) !! Per-draw mask for positive or negative infinity, excluding NaNs.
      logical, allocatable, intent(out) :: missing(:, :) !! Per-draw NaN mask; R's distinct NA payload is not represented.
      finite = ieee_is_finite(x)
      missing = ieee_is_nan(x)
      infinite = .not. finite .and. .not. missing
   end subroutine rvar_special_masks

   !> Append or replace the unnormalized log-weight variable of numeric draws.
   pure subroutine weight_draws_array(x, weights, values, log_scale, weight_variable, status)
      real(dp), intent(in) :: x(:, :, :) !! Numeric draws with dimensions (iteration, chain, variable).
      real(dp), intent(in) :: weights(:) !! One weight per pooled draw; no NaNs or positive infinity allowed.
      real(dp), allocatable, intent(out) :: values(:, :, :) !! Draws with an appended or replaced log-weight variable.
      logical, intent(in), optional :: log_scale !! Weights already on the log scale; default false.
      integer, intent(in), optional :: weight_variable !! Existing weight-column position; zero or absent appends a new column.
      integer, intent(out), optional :: status !! Zero on success, one for invalid weights, length, or column position.
      real(dp), allocatable :: lw(:)
      integer :: column, nv, i
      logical :: logarithmic
      allocate(values(0, 0, 0))
      if (present(status)) status = 1
      if (size(weights) /= size(x, 1) * size(x, 2)) return
      if (any(ieee_is_nan(weights))) return
      if (any(weights == ieee_value(0.0_dp, ieee_positive_inf))) return
      logarithmic = .false.
      if (present(log_scale)) logarithmic = log_scale
      column = 0
      if (present(weight_variable)) column = weight_variable
      if (column < 0 .or. column > size(x, 3)) return
      lw = weights
      if (.not. logarithmic) then
         if (any(weights < 0.0_dp)) return
         do i = 1, size(weights)
            if (weights(i) == 0.0_dp) then
               lw(i) = ieee_value(0.0_dp, ieee_negative_inf)
            else
               lw(i) = log(weights(i))
            end if
         end do
      end if
      nv = size(x, 3)
      if (column == 0) then
         nv = nv + 1
         column = nv
      end if
      deallocate(values)
      allocate(values(size(x, 1), size(x, 2), nv))
      values(:, :, :size(x, 3)) = x
      values(:, :, column) = reshape(lw, [size(x, 1), size(x, 2)])
      if (present(status)) status = 0
   end subroutine weight_draws_array

   !> Numerical weights.draws behavior, with explicit column metadata.
   pure subroutine extract_draws_weights(x, weight_variable, values, log_scale, normalize, status)
      real(dp), intent(in) :: x(:, :, :) !! Numeric draws containing an unnormalized log-weight column.
      integer, intent(in) :: weight_variable !! One-based column storing log weights.
      real(dp), allocatable, intent(out) :: values(:) !! Pooled normalized weights by default, empty on invalid column.
      logical, intent(in), optional :: log_scale !! Return log weights when true; default false returns linear weights.
      logical, intent(in), optional :: normalize !! Normalize by log-sum-exp when true (default).
      integer, intent(out), optional :: status !! Zero on success, one for an invalid column or stored NaN/positive infinity.
      logical :: logarithmic, normalized
      integer :: code
      call extract_draws_variable(x, weight_variable, values, code)
      if (present(status)) status = 1
      if (code /= 0) return
      if (any(ieee_is_nan(values)) .or. any(values == ieee_value(0.0_dp, ieee_positive_inf))) then
         values = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      logarithmic = .false.
      if (present(log_scale)) logarithmic = log_scale
      normalized = .true.
      if (present(normalize)) normalized = normalize
      if (normalized) values = values - r_log_sum_exp(values)
      if (.not. logarithmic) values = exp(values)
      if (present(status)) status = 0
   end subroutine extract_draws_weights

   pure subroutine rvar_over_draws_summary(x, method, values, na_rm, constant, status)
      real(dp), intent(in) :: x(:, :) !! Numeric rvar samples with dimensions (draw, element).
      character(len=*), intent(in) :: method !! Reduction across draws, using the same method names as rvar_summary.
      real(dp), allocatable, intent(out) :: values(:) !! One summary per rvar element.
      logical, intent(in), optional :: na_rm !! Remove NaNs independently per element; default false.
      real(dp), intent(in), optional :: constant !! Scale factor for MAD, default 1.4826.
      integer, intent(out), optional :: status !! Zero on success, one for an unknown reduction method.
      call rvar_summary(transpose(x), method, values, na_rm, constant, status)
   end subroutine rvar_over_draws_summary
end module posterior_draws
