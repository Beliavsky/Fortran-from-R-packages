! SPDX-License-Identifier: GPL-3.0-only
module neighbours
   use neighbours_kinds, only: dp, i8
   use neighbours_rng, only: rng_state, rng_seed, rng_uniform, rng_integer, &
      rng_shuffle, rng_sample
   implicit none
   private

   integer, parameter, public :: neighbours_ok = 0
   integer, parameter, public :: neighbours_invalid_input = 1
   integer, parameter, public :: neighbours_no_candidate = 2
   integer, parameter, public :: budget_zero_sum = 1
   integer, parameter, public :: budget_none = 2
   integer, parameter, public :: budget_range = 3

   type, public :: numeric_neighbour_config
      real(dp), allocatable :: lower(:), upper(:)
      integer, allocatable :: active(:)
      real(dp) :: stepsize = 0.0_dp
      logical :: random_step = .true.
      integer :: budget_mode = budget_zero_sum
      real(dp) :: budget_lower = -huge(1.0_dp)
      real(dp) :: budget_upper = huge(1.0_dp)
      logical :: update_ax = .false.
      real(dp), allocatable :: a(:, :)
   end type numeric_neighbour_config

   type, public :: logical_neighbour_config
      integer, allocatable :: active(:)
      integer :: stepsize = 1
      integer :: kmin = -1
      integer :: kmax = -1
   end type logical_neighbour_config

   public :: init_numeric_neighbour, init_logical_neighbour
   public :: numeric_neighbour, logical_neighbour
   public :: permute_neighbour, portfolio_5_10_40_neighbour
   public :: random_numeric_vector, random_numeric_vectors
   public :: random_logical_vector, random_logical_vectors
   public :: next_subset, compare_logical_vectors
   public :: rng_state, rng_seed, rng_uniform, rng_integer, rng_shuffle

   interface permute_neighbour
      module procedure permute_real
      module procedure permute_integer
      module procedure permute_logical
      module procedure permute_character
   end interface permute_neighbour
contains
   subroutine init_numeric_neighbour(config, n, stepsize, lower, upper, &
                                     random_step, budget_mode_in, budget, &
                                     active, update_ax, a, status)
      type(numeric_neighbour_config), intent(out) :: config
      integer, intent(in) :: n
      real(dp), intent(in) :: stepsize
      real(dp), intent(in), optional :: lower(:), upper(:), budget(:)
      logical, intent(in), optional :: random_step, update_ax
      integer, intent(in), optional :: budget_mode_in, active(:)
      real(dp), intent(in), optional :: a(:, :)
      integer, intent(out), optional :: status
      integer :: j, st

      st = neighbours_ok
      if (n <= 0 .or. stepsize < 0.0_dp) st = neighbours_invalid_input
      if (st /= neighbours_ok) then
         if (present(status)) status = st
         return
      end if

      allocate(config%lower(n), config%upper(n), config%active(n))
      config%lower = 0.0_dp
      config%upper = 1.0_dp
      config%active = [(j, j = 1, n)]
      config%stepsize = stepsize

      if (present(lower)) then
         if (size(lower) == 1) then
            config%lower = lower(1)
         else if (size(lower) == n) then
            config%lower = lower
         else
            st = neighbours_invalid_input
         end if
      end if
      if (present(upper)) then
         if (size(upper) == 1) then
            config%upper = upper(1)
         else if (size(upper) == n) then
            config%upper = upper
         else
            st = neighbours_invalid_input
         end if
      end if
      if (any(config%lower > config%upper)) st = neighbours_invalid_input

      if (present(random_step)) config%random_step = random_step
      if (present(budget_mode_in)) config%budget_mode = budget_mode_in
      if (config%budget_mode < budget_zero_sum .or. &
          config%budget_mode > budget_range) st = neighbours_invalid_input

      if (present(budget)) then
         if (size(budget) == 2) then
            config%budget_mode = budget_range
            config%budget_lower = budget(1)
            config%budget_upper = budget(2)
            if (budget(1) > budget(2)) st = neighbours_invalid_input
         else if (size(budget) == 1) then
            config%budget_mode = budget_zero_sum
         else
            st = neighbours_invalid_input
         end if
      end if

      if (present(active)) then
         if (size(active) < 1) then
            st = neighbours_invalid_input
         else if (any(active < 1) .or. any(active > n)) then
            st = neighbours_invalid_input
         else
            deallocate(config%active)
            allocate(config%active(size(active)))
            config%active = active
         end if
      end if

      if (present(update_ax)) config%update_ax = update_ax
      if (present(a)) then
         if (size(a, 2) /= n) then
            st = neighbours_invalid_input
         else
            allocate(config%a(size(a, 1), n))
            config%a = a
            config%update_ax = .true.
         end if
      else if (config%update_ax) then
         st = neighbours_invalid_input
      end if

      if (present(status)) status = st
   end subroutine init_numeric_neighbour

   subroutine init_logical_neighbour(config, n, stepsize, kmin, kmax, &
                                     active, status)
      type(logical_neighbour_config), intent(out) :: config
      integer, intent(in) :: n
      integer, intent(in), optional :: stepsize, kmin, kmax, active(:)
      integer, intent(out), optional :: status
      integer :: j, st

      st = neighbours_ok
      if (n <= 0) st = neighbours_invalid_input
      allocate(config%active(max(n, 0)))
      if (n > 0) config%active = [(j, j = 1, n)]
      if (present(stepsize)) config%stepsize = stepsize
      if (config%stepsize < 1) st = neighbours_invalid_input
      if (present(kmin)) config%kmin = kmin
      if (present(kmax)) config%kmax = kmax

      if ((config%kmin < 0) .neqv. (config%kmax < 0)) then
         st = neighbours_invalid_input
      else if (config%kmin >= 0) then
         if (config%kmin > config%kmax .or. config%kmax > n) &
            st = neighbours_invalid_input
      end if

      if (present(active)) then
         if (size(active) < 1 .or. any(active < 1) .or. any(active > n)) then
            st = neighbours_invalid_input
         else
            deallocate(config%active)
            allocate(config%active(size(active)))
            config%active = active
         end if
      end if
      if (present(status)) status = st
   end subroutine init_logical_neighbour

   subroutine numeric_neighbour(config, x, xn, rng, ax, status)
      type(numeric_neighbour_config), intent(in) :: config
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: xn(:)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(inout), optional :: ax(:)
      integer, intent(out), optional :: status
      integer, allocatable :: decrease(:), increase(:)
      integer :: i, j, nc, st
      real(dp) :: step, old_sum

      st = neighbours_ok
      xn = x
      if (.not. numeric_config_matches(config, size(x))) then
         st = neighbours_invalid_input
         if (present(status)) status = st
         return
      end if

      if (config%random_step) then
         step = config%stepsize * rng_uniform(rng)
      else
         step = config%stepsize
      end if

      if (config%budget_mode == budget_zero_sum) then
         call eligible_indices(x, config%lower, config%active, .true., decrease)
         call eligible_indices(x, config%upper, config%active, .false., increase)
         if (size(decrease) == 0 .or. size(increase) == 0) then
            st = neighbours_no_candidate
            if (present(status)) status = st
            return
         end if
         i = decrease(rng_integer(rng, size(decrease)))
         j = increase(rng_integer(rng, size(increase)))
         step = min(x(i) - config%lower(i), config%upper(j) - x(j), step)
         xn(i) = xn(i) - step
         xn(j) = xn(j) + step
      else
         old_sum = sum(x)
         if (rng_uniform(rng) < 0.5_dp) then
            step = -step
            call eligible_indices(x, config%lower, config%active, .true., decrease)
            nc = size(decrease)
            if (nc == 0) then
               st = neighbours_no_candidate
               if (present(status)) status = st
               return
            end if
            i = decrease(rng_integer(rng, nc))
            step = max(config%lower(i) - x(i), step)
            if (config%budget_mode == budget_range) &
               step = max(step, config%budget_lower - old_sum)
         else
            call eligible_indices(x, config%upper, config%active, .false., increase)
            nc = size(increase)
            if (nc == 0) then
               st = neighbours_no_candidate
               if (present(status)) status = st
               return
            end if
            i = increase(rng_integer(rng, nc))
            step = min(config%upper(i) - x(i), step)
            if (config%budget_mode == budget_range) &
               step = min(step, config%budget_upper - old_sum)
         end if
         xn(i) = xn(i) + step
      end if

      if (config%update_ax .and. present(ax)) then
         if (.not. allocated(config%a)) then
            st = neighbours_invalid_input
         else if (size(ax) /= size(config%a, 1)) then
            st = neighbours_invalid_input
         else
            ax = ax + matmul(config%a, xn - x)
         end if
      end if
      if (present(status)) status = st
   end subroutine numeric_neighbour

   logical function numeric_config_matches(config, n) result(ok)
      type(numeric_neighbour_config), intent(in) :: config
      integer, intent(in) :: n

      ok = allocated(config%lower) .and. allocated(config%upper) .and. &
           allocated(config%active)
      if (.not. ok) return
      ok = size(config%lower) == n .and. size(config%upper) == n
   end function numeric_config_matches

   subroutine eligible_indices(x, bound, active, for_decrease, indices)
      real(dp), intent(in) :: x(:), bound(:)
      integer, intent(in) :: active(:)
      logical, intent(in) :: for_decrease
      integer, allocatable, intent(out) :: indices(:)
      integer :: q, count

      count = 0
      do q = 1, size(active)
         if (for_decrease) then
            if (x(active(q)) > bound(active(q))) count = count + 1
         else
            if (x(active(q)) < bound(active(q))) count = count + 1
         end if
      end do
      allocate(indices(count))
      count = 0
      do q = 1, size(active)
         if (for_decrease) then
            if (x(active(q)) <= bound(active(q))) cycle
         else
            if (x(active(q)) >= bound(active(q))) cycle
         end if
         count = count + 1
         indices(count) = active(q)
      end do
   end subroutine eligible_indices

   subroutine logical_neighbour(config, x, xn, rng, status)
      type(logical_neighbour_config), intent(in) :: config
      logical, intent(in) :: x(:)
      logical, intent(out) :: xn(:)
      type(rng_state), intent(inout) :: rng
      integer, intent(out), optional :: status
      integer, allocatable :: pool(:), chosen(:), true_idx(:), false_idx(:)
      logical, allocatable :: xx(:)
      logical :: ok
      integer :: q, ntrue, st

      st = neighbours_ok
      xn = x
      if (.not. allocated(config%active) .or. config%stepsize < 1) then
         st = neighbours_invalid_input
         if (present(status)) status = st
         return
      end if
      if (any(config%active < 1) .or. any(config%active > size(x))) then
         st = neighbours_invalid_input
         if (present(status)) status = st
         return
      end if

      if (config%kmin < 0 .and. config%kmax < 0) then
         pool = config%active
         if (config%stepsize > size(pool)) then
            st = neighbours_invalid_input
         else
            allocate(chosen(config%stepsize))
            call rng_sample(rng, pool, config%stepsize, chosen, ok)
            if (.not. ok) then
               st = neighbours_invalid_input
            else
               do q = 1, size(chosen)
                  xn(chosen(q)) = .not. xn(chosen(q))
               end do
            end if
         end if
      else if (config%kmin == config%kmax) then
         allocate(xx(size(config%active)))
         xx = x(config%active)
         call logical_positions(xx, .true., true_idx)
         call logical_positions(xx, .false., false_idx)
         if (config%stepsize > size(true_idx) .or. &
             config%stepsize > size(false_idx)) then
            st = neighbours_no_candidate
         else
            allocate(chosen(config%stepsize))
            call rng_sample(rng, true_idx, config%stepsize, chosen, ok)
            do q = 1, size(chosen)
               xx(chosen(q)) = .false.
            end do
            call rng_sample(rng, false_idx, config%stepsize, chosen, ok)
            do q = 1, size(chosen)
               xx(chosen(q)) = .true.
            end do
            xn(config%active) = xx
         end if
      else
         ! The upstream kmin < kmax branch operates on the full vector and
         ! ignores active. Preserve that source behavior.
         call logical_positions(x, .true., true_idx)
         call logical_positions(x, .false., false_idx)
         ntrue = size(true_idx)
         allocate(chosen(config%stepsize))
         if (ntrue == config%kmax) then
            if (config%stepsize > size(true_idx)) then
               st = neighbours_no_candidate
            else
               call rng_sample(rng, true_idx, config%stepsize, chosen, ok)
               xn(chosen) = .false.
            end if
         else if (ntrue > config%kmin) then
            pool = [(q, q = 1, size(x))]
            if (config%stepsize > size(pool)) then
               st = neighbours_invalid_input
            else
               call rng_sample(rng, pool, config%stepsize, chosen, ok)
               do q = 1, size(chosen)
                  xn(chosen(q)) = .not. xn(chosen(q))
               end do
            end if
         else
            if (config%stepsize > size(false_idx)) then
               st = neighbours_no_candidate
            else
               call rng_sample(rng, false_idx, config%stepsize, chosen, ok)
               xn(chosen) = .true.
            end if
         end if
      end if
      if (present(status)) status = st
   end subroutine logical_neighbour

   subroutine logical_positions(x, value, indices)
      logical, intent(in) :: x(:), value
      integer, allocatable, intent(out) :: indices(:)
      integer :: j, ncount

      ncount = count(x .eqv. value)
      allocate(indices(ncount))
      ncount = 0
      do j = 1, size(x)
         if (x(j) .eqv. value) then
            ncount = ncount + 1
            indices(ncount) = j
         end if
      end do
   end subroutine logical_positions

   subroutine permutation_indices(rng, n, stepsize, indices, status)
      type(rng_state), intent(inout) :: rng
      integer, intent(in) :: n
      integer, intent(in), optional :: stepsize
      integer, allocatable, intent(out) :: indices(:)
      integer, intent(out) :: status
      integer, allocatable :: pool(:)
      integer :: j, k
      logical :: ok

      k = 2
      if (present(stepsize)) then
         if (stepsize /= 1) k = stepsize
      end if
      if (k < 2 .or. k > n) then
         status = neighbours_invalid_input
         allocate(indices(0))
         return
      end if
      pool = [(j, j = 1, n)]
      allocate(indices(k))
      call rng_sample(rng, pool, k, indices, ok)
      if (ok) then
         status = neighbours_ok
      else
         status = neighbours_invalid_input
      end if
   end subroutine permutation_indices

   subroutine permute_real(x, rng, stepsize, status)
      real(dp), intent(inout) :: x(:)
      type(rng_state), intent(inout) :: rng
      integer, intent(in), optional :: stepsize
      integer, intent(out), optional :: status
      integer, allocatable :: idx(:), order(:)
      real(dp), allocatable :: tmp(:)
      integer :: j, st

      call permutation_indices(rng, size(x), stepsize, idx, st)
      if (st == neighbours_ok) then
         if (size(idx) == 2) then
            tmp = [x(idx(1)), x(idx(2))]
            x(idx(1)) = tmp(2)
            x(idx(2)) = tmp(1)
         else
            order = [(j, j = 1, size(idx))]
            call rng_shuffle(rng, order)
            allocate(tmp(size(idx)))
            tmp = x(idx)
            x(idx) = tmp(order)
         end if
      end if
      if (present(status)) status = st
   end subroutine permute_real

   subroutine permute_integer(x, rng, stepsize, status)
      integer, intent(inout) :: x(:)
      type(rng_state), intent(inout) :: rng
      integer, intent(in), optional :: stepsize
      integer, intent(out), optional :: status
      integer, allocatable :: idx(:), order(:), tmp(:)
      integer :: j, st

      call permutation_indices(rng, size(x), stepsize, idx, st)
      if (st == neighbours_ok) then
         if (size(idx) == 2) then
            allocate(tmp(2))
            tmp = x(idx)
            x(idx(1)) = tmp(2)
            x(idx(2)) = tmp(1)
         else
            order = [(j, j = 1, size(idx))]
            call rng_shuffle(rng, order)
            allocate(tmp(size(idx)))
            tmp = x(idx)
            x(idx) = tmp(order)
         end if
      end if
      if (present(status)) status = st
   end subroutine permute_integer

   subroutine permute_logical(x, rng, stepsize, status)
      logical, intent(inout) :: x(:)
      type(rng_state), intent(inout) :: rng
      integer, intent(in), optional :: stepsize
      integer, intent(out), optional :: status
      integer, allocatable :: idx(:), order(:)
      logical, allocatable :: tmp(:)
      integer :: j, st

      call permutation_indices(rng, size(x), stepsize, idx, st)
      if (st == neighbours_ok) then
         order = [(j, j = 1, size(idx))]
         if (size(idx) == 2) then
            order = [2, 1]
         else
            call rng_shuffle(rng, order)
         end if
         allocate(tmp(size(idx)))
         tmp = x(idx)
         x(idx) = tmp(order)
      end if
      if (present(status)) status = st
   end subroutine permute_logical

   subroutine permute_character(x, rng, stepsize, status)
      character(len=*), intent(inout) :: x(:)
      type(rng_state), intent(inout) :: rng
      integer, intent(in), optional :: stepsize
      integer, intent(out), optional :: status
      integer, allocatable :: idx(:), order(:)
      character(len=len(x)), allocatable :: tmp(:)
      integer :: j, st

      call permutation_indices(rng, size(x), stepsize, idx, st)
      if (st == neighbours_ok) then
         order = [(j, j = 1, size(idx))]
         if (size(idx) == 2) then
            order = [2, 1]
         else
            call rng_shuffle(rng, order)
         end if
         allocate(tmp(size(idx)))
         tmp = x(idx)
         x(idx) = tmp(order)
      end if
      if (present(status)) status = st
   end subroutine permute_character

   subroutine portfolio_5_10_40_neighbour(x, xn, rng, kmax, status)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: xn(:)
      type(rng_state), intent(inout) :: rng
      integer, intent(in), optional :: kmax
      integer, intent(out), optional :: status
      integer, allocatable :: to_sell(:), to_buy(:)
      integer :: i, j, k, km, st
      real(dp), parameter :: wmax = 0.05_dp, wmax2 = 0.10_dp
      real(dp), parameter :: max_sum_l = 0.40_dp
      real(dp) :: eps, sum_l

      st = neighbours_ok
      xn = x
      km = 33
      if (present(kmax)) km = kmax
      k = count(abs(x) > 0.0_dp)
      eps = rng_uniform(rng) * 0.005_dp

      call positions_positive(x, to_sell)
      if (k == km) then
         call positions_between(x, 0.0_dp, wmax2, .true., to_buy)
      else
         call positions_less(x, wmax2, to_buy)
      end if
      if (size(to_sell) == 0 .or. size(to_buy) == 0) then
         st = neighbours_no_candidate
         if (present(status)) status = st
         return
      end if
      sum_l = sum(pack(x, x > wmax))
      i = to_sell(rng_integer(rng, size(to_sell)))
      j = to_buy(rng_integer(rng, size(to_buy)))

      if (x(j) < wmax) then
         eps = min(eps, wmax - x(j), x(i))
      else if (x(j) > wmax) then
         eps = min(eps, wmax2 - x(j), x(i), max(0.0_dp, max_sum_l - sum_l))
      else
         eps = min(eps, wmax2 - x(j), x(i), &
                   max(0.0_dp, max_sum_l - sum_l - x(j)))
      end if
      xn(i) = xn(i) - eps
      xn(j) = xn(j) + eps
      if (present(status)) status = st
   end subroutine portfolio_5_10_40_neighbour

   subroutine positions_positive(x, indices)
      real(dp), intent(in) :: x(:)
      integer, allocatable, intent(out) :: indices(:)
      integer :: j, q
      allocate(indices(count(x > 0.0_dp)))
      q = 0
      do j = 1, size(x)
         if (x(j) > 0.0_dp) then
            q = q + 1
            indices(q) = j
         end if
      end do
   end subroutine positions_positive

   subroutine positions_less(x, upper, indices)
      real(dp), intent(in) :: x(:), upper
      integer, allocatable, intent(out) :: indices(:)
      integer :: j, q
      allocate(indices(count(x < upper)))
      q = 0
      do j = 1, size(x)
         if (x(j) < upper) then
            q = q + 1
            indices(q) = j
         end if
      end do
   end subroutine positions_less

   subroutine positions_between(x, lower, upper, lower_open, indices)
      real(dp), intent(in) :: x(:), lower, upper
      logical, intent(in) :: lower_open
      integer, allocatable, intent(out) :: indices(:)
      logical, allocatable :: mask(:)
      integer :: j, q

      allocate(mask(size(x)))
      if (lower_open) then
         mask = x > lower .and. x < upper
      else
         mask = x >= lower .and. x < upper
      end if
      allocate(indices(count(mask)))
      q = 0
      do j = 1, size(x)
         if (mask(j)) then
            q = q + 1
            indices(q) = j
         end if
      end do
   end subroutine positions_between

   subroutine random_logical_vector(x, rng, kmin, kmax, status)
      logical, intent(out) :: x(:)
      type(rng_state), intent(inout) :: rng
      integer, intent(in), optional :: kmin, kmax
      integer, intent(out), optional :: status
      integer, allocatable :: pool(:), chosen(:)
      integer :: lo, hi, k, j, st
      logical :: ok

      lo = 0
      hi = size(x)
      if (present(kmin)) lo = kmin
      if (present(kmax)) hi = kmax
      st = neighbours_ok
      x = .false.
      if (lo < 0 .or. hi > size(x) .or. lo > hi) then
         st = neighbours_invalid_input
      else
         if (lo == hi) then
            k = lo
         else
            k = lo + rng_integer(rng, hi - lo + 1) - 1
         end if
         pool = [(j, j = 1, size(x))]
         allocate(chosen(k))
         call rng_sample(rng, pool, k, chosen, ok)
         if (k > 0) x(chosen) = .true.
      end if
      if (present(status)) status = st
   end subroutine random_logical_vector

   subroutine random_logical_vectors(x, rng, kmin, kmax, status)
      logical, intent(out) :: x(:, :)
      type(rng_state), intent(inout) :: rng
      integer, intent(in), optional :: kmin, kmax
      integer, intent(out), optional :: status
      integer :: j, st, local_status

      st = neighbours_ok
      do j = 1, size(x, 2)
         call random_logical_vector(x(:, j), rng, kmin, kmax, local_status)
         if (local_status /= neighbours_ok) st = local_status
      end do
      if (present(status)) status = st
   end subroutine random_logical_vectors

   subroutine random_numeric_vector(x, rng, lower, upper, kmin, kmax, status)
      real(dp), intent(out) :: x(:)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in), optional :: lower, upper
      integer, intent(in), optional :: kmin, kmax
      integer, intent(out), optional :: status
      integer, allocatable :: pool(:), chosen(:)
      real(dp) :: lo_value, hi_value
      integer :: lo, hi, nzero, j, st
      logical :: ok

      lo_value = 0.0_dp
      hi_value = 1.0_dp
      if (present(lower)) lo_value = lower
      if (present(upper)) hi_value = upper
      st = neighbours_ok
      if (lo_value > hi_value) st = neighbours_invalid_input
      do j = 1, size(x)
         x(j) = lo_value + (hi_value - lo_value) * rng_uniform(rng)
      end do
      if (present(kmin) .or. present(kmax)) then
         lo = 0
         hi = size(x)
         if (present(kmin)) lo = kmin
         if (present(kmax)) hi = kmax
         if (lo < 0 .or. hi > size(x) .or. lo > hi) then
            st = neighbours_invalid_input
         else
            if (lo == hi) then
               nzero = size(x) - lo
            else
               nzero = size(x) - (lo + rng_integer(rng, hi - lo + 1) - 1)
            end if
            pool = [(j, j = 1, size(x))]
            allocate(chosen(nzero))
            call rng_sample(rng, pool, nzero, chosen, ok)
            if (nzero > 0) x(chosen) = 0.0_dp
         end if
      end if
      if (present(status)) status = st
   end subroutine random_numeric_vector

   subroutine random_numeric_vectors(x, rng, lower, upper, kmin, kmax, status)
      real(dp), intent(out) :: x(:, :)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in), optional :: lower, upper
      integer, intent(in), optional :: kmin, kmax
      integer, intent(out), optional :: status
      integer :: j, st, local_status

      st = neighbours_ok
      do j = 1, size(x, 2)
         call random_numeric_vector(x(:, j), rng, lower, upper, kmin, kmax, &
                                    local_status)
         if (local_status /= neighbours_ok) st = local_status
      end do
      if (present(status)) status = st
   end subroutine random_numeric_vectors

   subroutine next_subset(a, n, k, has_next, status)
      integer, intent(inout) :: a(:)
      integer, intent(in) :: n, k
      logical, intent(out) :: has_next
      integer, intent(out), optional :: status
      integer :: h, h1, k1, m1, st

      st = neighbours_ok
      has_next = .false.
      if (k < 1 .or. n < k .or. size(a) /= k) then
         st = neighbours_invalid_input
         if (present(status)) status = st
         return
      end if
      if (a(1) == n - k + 1) then
         if (present(status)) status = st
         return
      end if
      h = 0
      do k1 = 1, k
         if (a(k + 1 - k1) /= n + 1 - k1) then
            h = k1
            exit
         end if
      end do
      if (h == 0) then
         if (present(status)) status = neighbours_invalid_input
         return
      end if
      m1 = a(k + 1 - h)
      do h1 = 1, h
         a(k + h1 - h) = m1 + h1
      end do
      has_next = .true.
      if (present(status)) status = st
   end subroutine next_subset

   subroutine compare_logical_vectors(vectors, distances, difference_mask, status)
      logical, intent(in) :: vectors(:, :)
      integer, allocatable, intent(out) :: distances(:)
      logical, allocatable, intent(out), optional :: difference_mask(:, :)
      integer, intent(out), optional :: status
      integer :: j

      if (size(vectors, 2) < 1) then
         allocate(distances(0))
         if (present(difference_mask)) allocate(difference_mask(size(vectors, 1), 0))
         if (present(status)) status = neighbours_invalid_input
         return
      end if
      allocate(distances(max(0, size(vectors, 2) - 1)))
      if (present(difference_mask)) &
         allocate(difference_mask(size(vectors, 1), max(0, size(vectors, 2) - 1)))
      do j = 2, size(vectors, 2)
         distances(j - 1) = count(vectors(:, j - 1) .neqv. vectors(:, j))
         if (present(difference_mask)) &
            difference_mask(:, j - 1) = vectors(:, j - 1) .neqv. vectors(:, j)
      end do
      if (present(status)) status = neighbours_ok
   end subroutine compare_logical_vectors
end module neighbours
