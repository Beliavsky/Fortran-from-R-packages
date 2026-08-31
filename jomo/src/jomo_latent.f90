! Latent-normal categorical augmentation and conditional imputation for jomo.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_latent
   use jomo_kinds, only : dp
   use jomo_rng, only : rng_state
   use jomo_distributions, only : mvnormal_sample
   use jomo_linalg, only : solve_spd, symmetrize
   implicit none
   private

   public :: latent_dimension
   public :: initialize_latent_data
   public :: set_categorical_covariance
   public :: decode_categories
   public :: update_observed_categories
   public :: impute_missing_rows

contains

   pure integer function latent_dimension(n_con, n_levels) result(p)
      integer, intent(in) :: n_con !! Number of continuous variables represented directly in the latent response vector.
      integer, intent(in) :: n_levels(:) !! Numbers of observed levels for categorical variables; each value must be at least two.

      if (n_con < 0) error stop "latent_dimension: n_con must be nonnegative"
      if (any(n_levels < 2)) error stop "latent_dimension: categorical variables need at least two levels"
      p = n_con + sum(n_levels - 1)
   end function latent_dimension

   subroutine initialize_latent_data(y_con, con_observed, y_cat, cat_observed, n_levels, latent, observed)
      real(dp), intent(in) :: y_con(:, :) !! Continuous responses, shape n by n_con; values flagged missing are ignored.
      logical, intent(in) :: con_observed(:, :) !! True where the corresponding continuous response is observed.
      integer, intent(in) :: y_cat(:, :) !! Categorical responses 1..K, shape n by n_cat; masked values are ignored.
      logical, intent(in) :: cat_observed(:, :) !! True where the corresponding categorical response is observed.
      integer, intent(in) :: n_levels(:) !! Number of levels K for each categorical response.
      real(dp), allocatable, intent(out) :: latent(:, :) !! Initialized expanded latent response matrix, shape n by p.
      logical, allocatable, intent(out) :: observed(:, :) !! Expanded mask for observed latent/continuous components.
      integer :: n
      integer :: n_con
      integer :: n_cat
      integer :: p
      integer :: i
      integer :: j
      integer :: k
      integer :: pos
      real(dp) :: mean_value

      n = size(y_con, 1)
      if (size(y_cat, 1) /= n) error stop "initialize_latent_data: row mismatch"
      if (any(shape(con_observed) /= shape(y_con))) error stop "initialize_latent_data: continuous mask shape mismatch"
      if (any(shape(cat_observed) /= shape(y_cat))) error stop "initialize_latent_data: categorical mask shape mismatch"
      n_con = size(y_con, 2)
      n_cat = size(y_cat, 2)
      if (size(n_levels) /= n_cat) error stop "initialize_latent_data: n_levels size mismatch"
      if (any(n_levels < 2)) error stop "initialize_latent_data: categorical variables need at least two levels"
      p = latent_dimension(n_con, n_levels)
      allocate(latent(n, p), observed(n, p))
      latent = 0.0_dp
      observed = .false.

      do j = 1, n_con
         if (count(con_observed(:, j)) > 0) then
            mean_value = sum(y_con(:, j), mask=con_observed(:, j)) / real(count(con_observed(:, j)), dp)
         else
            mean_value = 0.0_dp
         end if
         do i = 1, n
            if (con_observed(i, j)) then
               latent(i, j) = y_con(i, j)
               observed(i, j) = .true.
            else
               latent(i, j) = mean_value
            end if
         end do
      end do

      pos = n_con + 1
      do j = 1, n_cat
         do i = 1, n
            if (cat_observed(i, j)) then
               if (y_cat(i, j) < 1 .or. y_cat(i, j) > n_levels(j)) &
                  error stop "initialize_latent_data: category outside 1..K"
               observed(i, pos:pos + n_levels(j) - 2) = .true.
               latent(i, pos:pos + n_levels(j) - 2) = -0.5_dp
               if (y_cat(i, j) < n_levels(j)) then
                  k = pos + y_cat(i, j) - 1
                  latent(i, k) = 0.5_dp
               end if
            end if
         end do
         pos = pos + n_levels(j) - 1
      end do
   end subroutine initialize_latent_data

   pure subroutine set_categorical_covariance(omega, n_con, n_levels, fixed)
      real(dp), intent(inout) :: omega(:, :) !! Latent covariance; category blocks receive jomo identification values.
      integer, intent(in) :: n_con !! Number of leading continuous latent components.
      integer, intent(in) :: n_levels(:) !! Numbers of levels for categorical responses.
      logical, intent(out) :: fixed(:, :) !! Mask marking covariance elements fixed by categorical identification constraints.
      integer :: pos
      integer :: j
      integer :: k
      integer :: m

      if (size(omega, 1) /= size(omega, 2)) error stop "set_categorical_covariance: omega must be square"
      if (any(shape(fixed) /= shape(omega))) error stop "set_categorical_covariance: fixed mask shape mismatch"
      fixed = .false.
      pos = n_con + 1
      do m = 1, size(n_levels)
         do j = 0, n_levels(m) - 2
            do k = 0, n_levels(m) - 2
               if (j == k) then
                  omega(pos + j, pos + k) = 1.0_dp
               else
                  omega(pos + j, pos + k) = 0.5_dp
               end if
               fixed(pos + j, pos + k) = .true.
            end do
         end do
         pos = pos + n_levels(m) - 1
      end do
   end subroutine set_categorical_covariance

   pure subroutine decode_categories(latent, n_con, n_levels, categories)
      real(dp), intent(in) :: latent(:, :) !! Expanded latent response matrix, shape n by p.
      integer, intent(in) :: n_con !! Number of leading continuous latent components.
      integer, intent(in) :: n_levels(:) !! Number of levels for each categorical response.
      integer, intent(out) :: categories(:, :) !! Decoded n by n_cat categories; final level is the zero-reference class.
      integer :: i
      integer :: j
      integer :: k
      integer :: pos
      integer :: imax
      real(dp) :: vmax

      if (size(categories, 1) /= size(latent, 1) .or. size(categories, 2) /= size(n_levels)) &
         error stop "decode_categories: output shape mismatch"
      pos = n_con + 1
      do j = 1, size(n_levels)
         do i = 1, size(latent, 1)
            imax = 1
            vmax = latent(i, pos)
            do k = 2, n_levels(j) - 1
               if (latent(i, pos + k - 1) > vmax) then
                  vmax = latent(i, pos + k - 1)
                  imax = k
               end if
            end do
            if (vmax > 0.0_dp) then
               categories(i, j) = imax
            else
               categories(i, j) = n_levels(j)
            end if
         end do
         pos = pos + n_levels(j) - 1
      end do
   end subroutine decode_categories

   subroutine update_observed_categories(rng, categories, cat_observed, n_levels, n_con, mean, omega, latent, failures)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for constrained latent-normal draws.
      integer, intent(in) :: categories(:, :) !! Observed categorical labels encoded 1..K, shape n by n_cat.
      logical, intent(in) :: cat_observed(:, :) !! True where category labels are observed and must constrain the latent draw.
      integer, intent(in) :: n_levels(:) !! Number of levels K for each categorical variable.
      integer, intent(in) :: n_con !! Number of leading continuous components in the latent vector.
      real(dp), intent(in) :: mean(:, :) !! Current row-specific latent means, shape n by p.
      real(dp), intent(in) :: omega(:, :) !! Current common latent covariance matrix, p by p.
      real(dp), intent(inout) :: latent(:, :) !! Latent response matrix updated in-place for observed categorical blocks.
      integer, intent(out) :: failures !! Number of categorical blocks with no accepted draw in 10000 attempts.
      integer :: n
      integer :: p
      integer :: pos
      integer :: m
      integer :: i
      integer :: j
      integer :: o
      integer :: info
      integer :: tries
      integer :: selected
      integer, allocatable :: other(:)
      real(dp), allocatable :: soo(:, :)
      real(dp), allocatable :: sbo(:, :)
      real(dp), allocatable :: sb(:, :)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: diff(:)
      real(dp), allocatable :: mb(:)
      real(dp), allocatable :: draw(:)
      logical :: accepted

      n = size(latent, 1)
      p = size(latent, 2)
      if (any(shape(mean) /= shape(latent))) error stop "update_observed_categories: mean shape mismatch"
      if (size(omega, 1) /= p .or. size(omega, 2) /= p) error stop "update_observed_categories: omega shape mismatch"
      if (size(categories, 1) /= n .or. size(categories, 2) /= size(n_levels)) &
         error stop "update_observed_categories: category shape mismatch"
      if (any(shape(cat_observed) /= shape(categories))) error stop "update_observed_categories: mask shape mismatch"

      failures = 0
      pos = n_con + 1
      do m = 1, size(n_levels)
         j = n_levels(m) - 1
         o = p - j
         allocate(sb(j, j), mb(j), draw(j))
         sb = omega(pos:pos + j - 1, pos:pos + j - 1)
         if (o > 0) then
            allocate(other(o), soo(o, o), sbo(j, o), v(o, j), diff(o))
            call complement_indices(p, pos, pos + j - 1, other)
            soo = omega(other, other)
            sbo = omega(pos:pos + j - 1, other)
            call solve_spd(soo, transpose(sbo), v, info)
            if (info /= 0) error stop "update_observed_categories: conditional covariance solve failed"
            sb = sb - matmul(sbo, v)
            call symmetrize(sb)
         end if

         do i = 1, n
            if (.not. cat_observed(i, m)) cycle
            mb = mean(i, pos:pos + j - 1)
            if (o > 0) then
               diff = latent(i, other) - mean(i, other)
               mb = mb + matmul(transpose(v), diff)
            end if
            selected = categories(i, m)
            accepted = .false.
            do tries = 1, 10000
               call mvnormal_sample(rng, mb, sb, draw, info)
               if (info /= 0) exit
               if (minval(draw) <= -3.0_dp .or. maxval(draw) >= 4.0_dp) cycle
               if (selected == n_levels(m)) then
                  accepted = maxval(draw) < 0.0_dp
               else
                  accepted = maxloc(draw, dim=1) == selected .and. maxval(draw) > 0.0_dp
               end if
               if (accepted) exit
            end do
            if (accepted) then
               latent(i, pos:pos + j - 1) = draw
            else
               failures = failures + 1
            end if
         end do
         if (o > 0) deallocate(other, soo, sbo, v, diff)
         deallocate(sb, mb, draw)
         pos = pos + j
      end do
   end subroutine update_observed_categories

   subroutine impute_missing_rows(rng, observed, mean, covariance, values, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for conditional multivariate-normal imputations.
      logical, intent(in) :: observed(:, :) !! Fixed observation mask; false components are redrawn at every call.
      real(dp), intent(in) :: mean(:, :) !! Current row-specific multivariate-normal means, same shape as values.
      real(dp), intent(in) :: covariance(:, :) !! Common covariance matrix of the complete latent response vector.
      real(dp), intent(inout) :: values(:, :) !! Complete latent response matrix whose missing components are updated in-place.
      integer, intent(out) :: info !! Zero on success; positive if a conditional covariance factorization fails.
      integer :: i
      integer :: p
      integer :: nm
      integer :: no
      integer :: k
      integer :: solve_info
      integer, allocatable :: miss(:)
      integer, allocatable :: obs(:)
      real(dp), allocatable :: smm(:, :)
      real(dp), allocatable :: soo(:, :)
      real(dp), allocatable :: smo(:, :)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: conditional_mean(:)
      real(dp), allocatable :: conditional_cov(:, :)
      real(dp), allocatable :: diff(:)
      real(dp), allocatable :: draw(:)

      if (any(shape(observed) /= shape(values))) error stop "impute_missing_rows: observation mask shape mismatch"
      if (any(shape(mean) /= shape(values))) error stop "impute_missing_rows: mean shape mismatch"
      p = size(values, 2)
      if (size(covariance, 1) /= p .or. size(covariance, 2) /= p) error stop "impute_missing_rows: covariance shape mismatch"
      info = 0

      do i = 1, size(values, 1)
         nm = count(.not. observed(i, :))
         if (nm == 0) cycle
         no = p - nm
         allocate(miss(nm), smm(nm, nm), conditional_mean(nm), conditional_cov(nm, nm), draw(nm))
         k = 0
         call collect_indices(.not. observed(i, :), miss)
         smm = covariance(miss, miss)
         conditional_mean = mean(i, miss)
         conditional_cov = smm
         if (no > 0) then
            allocate(obs(no), soo(no, no), smo(nm, no), v(no, nm), diff(no))
            call collect_indices(observed(i, :), obs)
            soo = covariance(obs, obs)
            smo = covariance(miss, obs)
            call solve_spd(soo, transpose(smo), v, solve_info)
            if (solve_info /= 0) then
               info = solve_info
               return
            end if
            diff = values(i, obs) - mean(i, obs)
            conditional_mean = conditional_mean + matmul(transpose(v), diff)
            conditional_cov = smm - matmul(smo, v)
            call symmetrize(conditional_cov)
            deallocate(obs, soo, smo, v, diff)
         end if
         call mvnormal_sample(rng, conditional_mean, conditional_cov, draw, solve_info)
         if (solve_info /= 0) then
            info = solve_info
            return
         end if
         values(i, miss) = draw
         deallocate(miss, smm, conditional_mean, conditional_cov, draw)
      end do
   end subroutine impute_missing_rows

   pure subroutine complement_indices(p, first, last, idx)
      integer, intent(in) :: p !! Total number of indices in the source vector.
      integer, intent(in) :: first !! First one-based index excluded from the complement.
      integer, intent(in) :: last !! Last one-based index excluded from the complement.
      integer, intent(out) :: idx(:) !! One-based indices outside the inclusive first:last block.
      integer :: i
      integer :: k

      k = 0
      do i = 1, p
         if (i >= first .and. i <= last) cycle
         k = k + 1
         idx(k) = i
      end do
   end subroutine complement_indices

   pure subroutine collect_indices(mask, idx)
      logical, intent(in) :: mask(:) !! Logical mask identifying indices to collect.
      integer, intent(out) :: idx(:) !! One-based positions of true mask entries, in increasing order.
      integer :: i
      integer :: k

      k = 0
      do i = 1, size(mask)
         if (mask(i)) then
            k = k + 1
            idx(k) = i
         end if
      end do
   end subroutine collect_indices

end module jomo_latent
