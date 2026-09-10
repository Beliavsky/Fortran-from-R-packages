! SPDX-License-Identifier: GPL-2.0-only
module marss_model_ops
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use marss_kinds, only : dp
   use r_linalg, only : spectral_radius
   use marss_types, only : marss_model, marss_constraint_block, marss_constraints
   use marss_utils, only : identity_matrix, is_symmetric
   use marss_parameters, only : marss_time_shape_valid, marss_b_at, marss_q_at, marss_r_at, marss_v0_effective
   use marss_covariance, only : is_positive_semidefinite
   implicit none
   private
   public :: marss_inits
   public :: marss_kemcheck_structural
   public :: covariance_diagonal_fixed_zero
   public :: marss_model_valid

contains

   pure subroutine marss_inits(y, nstate, model, tinitx)
      real(dp), intent(in) :: y(:, :) !! Observation matrix with variables by time, NaNs mark missing values.
      integer, intent(in) :: nstate !! Number of latent states to initialize.
      type(marss_model), intent(out) :: model !! Initialized time-invariant MARSS model.
      integer, intent(in), optional :: tinitx !! Initial-state convention: zero for x(0), one for x(1).
      integer :: i
      integer :: nobs
      integer :: t

      nobs = size(y, 1)
      allocate(model%y(nobs, size(y, 2)))
      allocate(model%b(nstate, nstate), model%u(nstate), model%q(nstate, nstate))
      allocate(model%z(nobs, nstate), model%a(nobs), model%r(nobs, nobs))
      allocate(model%x0(nstate), model%v0(nstate, nstate))
      model%y = y
      model%b = identity_matrix(nstate)
      model%u = 0.0_dp
      model%q = 0.05_dp * identity_matrix(nstate)
      model%z = 0.0_dp
      do i = 1, min(nobs, nstate)
         model%z(i, i) = 1.0_dp
      end do
      if (nstate == 1 .and. nobs > 1) model%z(:, 1) = 1.0_dp
      model%a = 0.0_dp
      model%r = 0.05_dp * identity_matrix(nobs)
      model%x0 = 0.0_dp
      do i = 1, min(nobs, nstate)
         do t = 1, size(y, 2)
            if (.not. ieee_is_nan(y(i, t))) then
               model%x0(i) = y(i, t)
               exit
            end if
         end do
      end do
      model%v0 = 5.0_dp * identity_matrix(nstate)
      model%tinitx = 0
      if (present(tinitx)) model%tinitx = tinitx
   end subroutine marss_inits

   pure logical function marss_model_valid(model) result(ok)
      type(marss_model), intent(in) :: model !! Model whose dimensions, time slices, and covariance structure are checked.
      integer :: m
      integer :: n
      integer :: t
      integer :: tt

      ok = .false.
      if (.not. allocated(model%y)) return
      if (.not. allocated(model%b)) return
      if (.not. allocated(model%u)) return
      if (.not. allocated(model%q)) return
      if (.not. allocated(model%z)) return
      if (.not. allocated(model%a)) return
      if (.not. allocated(model%r)) return
      if (.not. allocated(model%x0)) return
      if (.not. allocated(model%v0)) return
      n = size(model%y, 1)
      tt = size(model%y, 2)
      m = size(model%b, 1)
      if (m < 1 .or. n < 1 .or. tt < 1) return
      if (size(model%b, 2) /= m) return
      if (size(model%u) /= m) return
      if (any(shape(model%q) /= [m, m])) return
      if (any(shape(model%z) /= [n, m])) return
      if (size(model%a) /= n) return
      if (any(shape(model%r) /= [n, n])) return
      if (size(model%x0) /= m) return
      if (any(shape(model%v0) /= [m, m])) return
      if (.not. marss_time_shape_valid(model, tt)) return
      if (allocated(model%q_noise) .and. allocated(model%q_t)) return
      if (allocated(model%r_noise) .and. allocated(model%r_t)) return
      if (.not. covariance_valid(model%q)) return
      if (.not. covariance_valid(model%r)) return
      if (.not. covariance_valid(model%v0)) return
      if (allocated(model%q_t)) then
         do t = 1, tt
            if (.not. covariance_valid(model%q_t(:, :, t))) return
         end do
      end if
      if (allocated(model%r_t)) then
         do t = 1, tt
            if (.not. covariance_valid(model%r_t(:, :, t))) return
         end do
      end if
      if (allocated(model%q_noise)) then
         if (.not. covariance_valid(model%q_noise)) return
         if (allocated(model%q_noise_t)) then
            do t = 1, tt
               if (.not. covariance_valid(model%q_noise_t(:, :, t))) return
            end do
         end if
      end if
      if (allocated(model%r_noise)) then
         if (.not. covariance_valid(model%r_noise)) return
         if (allocated(model%r_noise_t)) then
            do t = 1, tt
               if (.not. covariance_valid(model%r_noise_t(:, :, t))) return
            end do
         end if
      end if
      if (allocated(model%l) .neqv. allocated(model%v0_noise)) return
      if (allocated(model%l)) then
         if (size(model%l, 1) /= m) return
         if (size(model%v0_noise, 1) /= size(model%l, 2)) return
         if (size(model%v0_noise, 2) /= size(model%l, 2)) return
         if (.not. covariance_valid(model%v0_noise)) return
      end if
      if (.not. covariance_valid(marss_v0_effective(model))) return
      do t = 1, tt
         if (.not. covariance_valid(marss_q_at(model, t))) return
         if (.not. covariance_valid(marss_r_at(model, t))) return
      end do
      if (model%tinitx /= 0 .and. model%tinitx /= 1) return
      if (model%diffuse .and. model%tinitx /= 1) return
      ok = .true.
   end function marss_model_valid

   subroutine marss_kemcheck_structural(model, ok, info, constraints)
      type(marss_model), intent(in) :: model !! Numerical MARSS model checked for EM-compatible degeneracy restrictions.
      logical, intent(out) :: ok !! True when the model satisfies the implemented KEM validity conditions.
      integer, intent(out) :: info !! Zero on success; positive values identify the first failed KEM restriction.
      type(marss_constraints), intent(in), optional :: constraints !! Optional affine f+D*beta structure identifying free elements.
      real(dp), allocatable :: b_now(:, :)
      real(dp), allocatable :: q_now(:, :)
      real(dp), allocatable :: r_now(:, :)
      logical, allocatable :: b_adjacency(:, :)
      logical, allocatable :: b_adjacency_ref(:, :)
      logical, allocatable :: deterministic_now(:)
      logical, allocatable :: deterministic_ref(:)
      logical, allocatable :: indirectly_stochastic_now(:)
      logical, allocatable :: indirectly_stochastic_ref(:)
      logical, allocatable :: q_zero_ref(:)
      logical, allocatable :: q_zero(:)
      logical, allocatable :: r_zero(:)
      logical :: b_adjacency_changed
      logical :: check_b_adjacency
      logical :: deterministic_changed
      logical :: indirectly_stochastic_changed
      real(dp) :: radius
      real(dp) :: radius_tolerance
      integer :: eig_info
      integer :: i
      integer :: j
      integer :: m
      integer :: n
      integer :: t
      integer :: tt

      ok = .false.
      info = 1
      if (.not. marss_model_valid(model)) return
      tt = size(model%y, 2)
      if (tt <= 2) then
         info = 2
         return
      end if
      if (present(constraints)) then
         if (.not. kem_constraint_shapes_valid(model, constraints)) then
            info = 10
            return
         end if
      end if

      do t = 1, merge(tt, 1, allocated(model%b_t))
         if (present(constraints)) then
            if (matrix_slice_free(constraints%b, model%b, model%b_t, t)) cycle
         end if
         b_now = marss_b_at(model, t)
         call spectral_radius(b_now, radius, eig_info)
         if (eig_info /= 0) then
            info = 12
            return
         end if
         radius_tolerance = sqrt(epsilon(1.0_dp))
         if (radius > 1.0_dp + radius_tolerance) then
            info = 11
            return
         end if
      end do

      if (.not. present(constraints)) then
         ok = .true.
         info = 0
         return
      end if

      m = size(model%b, 1)
      n = size(model%z, 1)
      allocate(q_zero_ref(m), q_zero(m), r_zero(n))
      allocate(b_adjacency(m, m), b_adjacency_ref(m, m))
      allocate(deterministic_now(m), deterministic_ref(m))
      allocate(indirectly_stochastic_now(m), indirectly_stochastic_ref(m))
      b_adjacency_changed = .false.
      check_b_adjacency = .false.
      deterministic_changed = .false.
      indirectly_stochastic_changed = .false.

      do t = 1, tt
         q_now = marss_q_at(model, t)
         r_now = marss_r_at(model, t)
         do j = 1, m
            q_zero(j) = covariance_diagonal_fixed_zero(model, constraints, q_now(j, j), j, t, .true.)
         end do
         do i = 1, n
            r_zero(i) = covariance_diagonal_fixed_zero(model, constraints, r_now(i, i), i, t, .false.)
         end do

         if (t == 1) then
            q_zero_ref = q_zero
         else if (any(q_zero .neqv. q_zero_ref)) then
            info = 30
            return
         end if

         call potential_b_adjacency_at(model, constraints, t, b_adjacency)
         call kem_state_classes(b_adjacency, q_zero, deterministic_now, indirectly_stochastic_now)
         do j = 1, m
            if (.not. indirectly_stochastic_now(j)) cycle
            if (vector_element_free(constraints%u, model%u, model%u_t, j, t)) then
               info = 37
               return
            end if
            if (allocated(model%c_coef)) then
               if (coefficient_row_free(constraints%c, model%c_coef, model%c_coef_t, j, t)) then
                  info = 37
                  return
               end if
            end if
         end do
         if (t == 1) then
            b_adjacency_ref = b_adjacency
            deterministic_ref = deterministic_now
            indirectly_stochastic_ref = indirectly_stochastic_now
         else
            b_adjacency_changed = b_adjacency_changed .or. any(b_adjacency .neqv. b_adjacency_ref)
            deterministic_changed = deterministic_changed .or. any(deterministic_now .neqv. deterministic_ref)
            indirectly_stochastic_changed = indirectly_stochastic_changed .or. &
               any(indirectly_stochastic_now .neqv. indirectly_stochastic_ref)
         end if

         do i = 1, n
            if (.not. r_zero(i)) cycle
            if (matrix_row_free(constraints%z, model%z, model%z_t, i, t)) then
               info = 20
               return
            end if
            if (vector_element_free(constraints%a, model%a, model%a_t, i, t)) then
               info = 21
               return
            end if
            if (allocated(model%d_coef)) then
               if (coefficient_row_free(constraints%d, model%d_coef, model%d_coef_t, i, t)) then
                  info = 21
                  return
               end if
            end if
         end do

         do j = 1, m
            if (.not. q_zero(j)) cycle
            if (vector_element_free(constraints%u, model%u, model%u_t, j, t)) check_b_adjacency = .true.
            if (matrix_row_free(constraints%b, model%b, model%b_t, j, t)) then
               info = 31
               return
            end if
            do i = 1, n
               if (.not. r_zero(i)) cycle
               if (.not. potential_matrix_nonzero(constraints%z, model%z, model%z_t, i, j, t)) cycle
               if (vector_element_free(constraints%u, model%u, model%u_t, j, t)) then
                  info = 33
                  return
               end if
               if (allocated(model%c_coef)) then
                  if (coefficient_row_free(constraints%c, model%c_coef, model%c_coef_t, j, t)) then
                     info = 33
                     return
                  end if
               end if
            end do
         end do
      end do

      do j = 1, m
         if (.not. q_zero_ref(j)) cycle
         if (constraint_element_free(constraints%x0, j)) check_b_adjacency = .true.
      end do
      if (check_b_adjacency .and. b_adjacency_changed) then
         info = 34
         return
      end if
      if (deterministic_changed) then
         info = 35
         return
      end if
      if (indirectly_stochastic_changed) then
         info = 36
         return
      end if

      ok = .true.
      info = 0
   end subroutine marss_kemcheck_structural

   pure logical function kem_constraint_shapes_valid(model, constraints) result(ok)
      type(marss_model), intent(in) :: model !! Model defining flattened block sizes for the affine KEM checks.
      type(marss_constraints), intent(in) :: constraints !! Affine constraint blocks whose dimensions are validated.

      ok = constraint_block_shape_valid(constraints%b, active_matrix_size(model%b, model%b_t))
      if (.not. ok) return
      ok = constraint_block_shape_valid(constraints%u, active_vector_size(model%u, model%u_t))
      if (.not. ok) return
      ok = constraint_block_shape_valid(constraints%q, active_q_raw_size(model))
      if (.not. ok) return
      ok = constraint_block_shape_valid(constraints%z, active_matrix_size(model%z, model%z_t))
      if (.not. ok) return
      ok = constraint_block_shape_valid(constraints%a, active_vector_size(model%a, model%a_t))
      if (.not. ok) return
      ok = constraint_block_shape_valid(constraints%r, active_r_raw_size(model))
      if (.not. ok) return
      ok = constraint_block_shape_valid(constraints%x0, size(model%x0))
      if (.not. ok) return
      if (allocated(model%c_coef)) then
         ok = constraint_block_shape_valid(constraints%c, active_matrix_size(model%c_coef, model%c_coef_t))
         if (.not. ok) return
      else if (allocated(constraints%c%design)) then
         ok = .false.
         return
      end if
      if (allocated(model%d_coef)) then
         ok = constraint_block_shape_valid(constraints%d, active_matrix_size(model%d_coef, model%d_coef_t))
         if (.not. ok) return
      else if (allocated(constraints%d%design)) then
         ok = .false.
         return
      end if
      if (allocated(model%g)) then
         ok = constraint_block_shape_valid(constraints%g, active_matrix_size(model%g, model%g_t))
         if (.not. ok) return
      else if (allocated(constraints%g%design)) then
         ok = .false.
         return
      end if
      if (allocated(model%h)) then
         ok = constraint_block_shape_valid(constraints%h, active_matrix_size(model%h, model%h_t))
         if (.not. ok) return
      else if (allocated(constraints%h%design)) then
         ok = .false.
      end if
   end function kem_constraint_shapes_valid

   pure logical function constraint_block_shape_valid(block, nelem) result(ok)
      type(marss_constraint_block), intent(in) :: block !! Affine block checked for fixed, design, and start dimensions.
      integer, intent(in) :: nelem !! Number of flattened elements in the associated numerical block.

      ok = .true.
      if (.not. allocated(block%design)) then
         if (allocated(block%fixed) .or. allocated(block%start)) ok = .false.
         return
      end if
      if (.not. allocated(block%fixed) .or. .not. allocated(block%start)) then
         ok = .false.
         return
      end if
      if (size(block%fixed) /= nelem) ok = .false.
      if (size(block%design, 1) /= nelem) ok = .false.
      if (size(block%design, 2) /= size(block%start)) ok = .false.
   end function constraint_block_shape_valid

   pure integer function active_matrix_size(base, timed) result(nelem)
      real(dp), intent(in) :: base(:, :) !! Static matrix whose size is used when no time-indexed representation is allocated.
      real(dp), allocatable, intent(in) :: timed(:, :, :) !! Optional time-indexed representation of the same parameter block.

      if (allocated(timed)) then
         nelem = size(timed)
      else
         nelem = size(base)
      end if
   end function active_matrix_size

   pure integer function active_vector_size(base, timed) result(nelem)
      real(dp), intent(in) :: base(:) !! Static vector whose size is used when no time-indexed representation is allocated.
      real(dp), allocatable, intent(in) :: timed(:, :) !! Optional time-indexed representation of the same parameter block.

      if (allocated(timed)) then
         nelem = size(timed)
      else
         nelem = size(base)
      end if
   end function active_vector_size

   pure integer function active_q_raw_size(model) result(nelem)
      type(marss_model), intent(in) :: model !! Model whose raw Q or Q-noise affine block size is requested.

      if (allocated(model%q_noise)) then
         if (allocated(model%q_noise_t)) then
            nelem = size(model%q_noise_t)
         else
            nelem = size(model%q_noise)
         end if
      else
         nelem = active_matrix_size(model%q, model%q_t)
      end if
   end function active_q_raw_size

   pure integer function active_r_raw_size(model) result(nelem)
      type(marss_model), intent(in) :: model !! Model whose raw R or R-noise affine block size is requested.

      if (allocated(model%r_noise)) then
         if (allocated(model%r_noise_t)) then
            nelem = size(model%r_noise_t)
         else
            nelem = size(model%r_noise)
         end if
      else
         nelem = active_matrix_size(model%r, model%r_t)
      end if
   end function active_r_raw_size

   pure logical function covariance_diagonal_fixed_zero(model, constraints, value, component, t, process) result(is_zero)
      type(marss_model), intent(in) :: model !! Model supplying direct or loading-factor covariance representations.
      type(marss_constraints), intent(in) :: constraints !! Affine constraints used to distinguish fixed from estimated variances.
      real(dp), intent(in) :: value !! Effective covariance diagonal value at the requested time and component.
      integer, intent(in) :: component !! One-based state or observation component whose variance is tested.
      integer, intent(in) :: t !! One-based time slice for a time-indexed covariance.
      logical, intent(in) :: process !! True for process covariance Q; false for observation covariance R.
      real(dp) :: scale
      real(dp) :: tolerance
      integer :: idx
      integer :: nrow
      integer :: tslice

      scale = max(1.0_dp, abs(value))
      tolerance = 128.0_dp * epsilon(1.0_dp) * scale
      is_zero = abs(value) <= tolerance
      if (.not. is_zero) return

      if (process) then
         if (allocated(model%g) .and. allocated(model%q_noise)) then
            is_zero = .not. constraint_block_has_free(constraints%g) .and. &
               .not. constraint_block_has_free(constraints%q)
            return
         end if
         nrow = size(model%q, 1)
         tslice = merge(t, 1, allocated(model%q_t))
         idx = component + (component - 1) * nrow + (tslice - 1) * nrow * nrow
         is_zero = .not. constraint_element_free(constraints%q, idx)
      else
         if (allocated(model%h) .and. allocated(model%r_noise)) then
            is_zero = .not. constraint_block_has_free(constraints%h) .and. &
               .not. constraint_block_has_free(constraints%r)
            return
         end if
         nrow = size(model%r, 1)
         tslice = merge(t, 1, allocated(model%r_t))
         idx = component + (component - 1) * nrow + (tslice - 1) * nrow * nrow
         is_zero = .not. constraint_element_free(constraints%r, idx)
      end if
   end function covariance_diagonal_fixed_zero

   pure logical function constraint_block_has_free(block) result(has_free)
      type(marss_constraint_block), intent(in) :: block !! Affine block tested for any nonzero design coefficient.

      has_free = .false.
      if (.not. allocated(block%design)) return
      if (size(block%design, 2) == 0) return
      has_free = any(abs(block%design) > 128.0_dp * epsilon(1.0_dp))
   end function constraint_block_has_free

   pure logical function matrix_slice_free(block, base, timed, t) result(is_free)
      type(marss_constraint_block), intent(in) :: block !! Affine constraint block for the matrix time slice being inspected.
      real(dp), intent(in) :: base(:, :) !! Static matrix defining the row and column dimensions.
      real(dp), allocatable, intent(in) :: timed(:, :, :) !! Optional time-indexed matrix representation.
      integer, intent(in) :: t !! One-based time index used for a time-varying block.
      integer :: col
      integer :: idx
      integer :: row
      integer :: tslice

      is_free = .false.
      if (.not. allocated(block%design)) return
      tslice = merge(t, 1, allocated(timed))
      do col = 1, size(base, 2)
         do row = 1, size(base, 1)
            idx = row + (col - 1) * size(base, 1) + (tslice - 1) * size(base)
            if (constraint_element_free(block, idx)) then
               is_free = .true.
               return
            end if
         end do
      end do
   end function matrix_slice_free

   pure logical function matrix_row_free(block, base, timed, row, t) result(is_free)
      type(marss_constraint_block), intent(in) :: block !! Affine constraint block for the matrix being inspected.
      real(dp), intent(in) :: base(:, :) !! Static matrix used when no time-indexed representation is active.
      real(dp), allocatable, intent(in) :: timed(:, :, :) !! Optional time-indexed matrix representation.
      integer, intent(in) :: row !! One-based matrix row whose free elements are tested.
      integer, intent(in) :: t !! One-based time index used for a time-varying block.
      integer :: col
      integer :: idx
      integer :: tslice

      is_free = .false.
      if (.not. allocated(block%design)) return
      tslice = merge(t, 1, allocated(timed))
      do col = 1, size(base, 2)
         idx = row + (col - 1) * size(base, 1) + (tslice - 1) * size(base)
         if (constraint_element_free(block, idx)) then
            is_free = .true.
            return
         end if
      end do
   end function matrix_row_free

   pure logical function vector_element_free(block, base, timed, row, t) result(is_free)
      type(marss_constraint_block), intent(in) :: block !! Affine constraint block for the vector being inspected.
      real(dp), intent(in) :: base(:) !! Static vector used when no time-indexed representation is active.
      real(dp), allocatable, intent(in) :: timed(:, :) !! Optional time-indexed vector representation.
      integer, intent(in) :: row !! One-based vector element whose free status is tested.
      integer, intent(in) :: t !! One-based time index used for a time-varying block.
      integer :: idx
      integer :: tslice

      is_free = .false.
      if (.not. allocated(block%design)) return
      tslice = merge(t, 1, allocated(timed))
      idx = row + (tslice - 1) * size(base)
      is_free = constraint_element_free(block, idx)
   end function vector_element_free

   pure logical function coefficient_row_free(block, base, timed, row, t) result(is_free)
      type(marss_constraint_block), intent(in) :: block !! Affine constraint block for C or D coefficients.
      real(dp), intent(in) :: base(:, :) !! Static coefficient matrix defining row and covariate dimensions.
      real(dp), allocatable, intent(in) :: timed(:, :, :) !! Optional time-indexed coefficient representation.
      integer, intent(in) :: row !! One-based coefficient row associated with a state or observation component.
      integer, intent(in) :: t !! One-based time index used for a time-varying coefficient block.

      is_free = matrix_row_free(block, base, timed, row, t)
   end function coefficient_row_free

   pure subroutine potential_b_adjacency_at(model, constraints, t, adjacency)
      type(marss_model), intent(in) :: model !! Model supplying static or time-indexed B entries.
      type(marss_constraints), intent(in) :: constraints !! Affine B constraints used to identify potentially nonzero entries.
      integer, intent(in) :: t !! One-based time index at which the B adjacency pattern is constructed.
      logical, intent(out) :: adjacency(:, :) !! Potential B support matrix; true entries may be nonzero under the constraints.
      integer :: col
      integer :: row

      adjacency = .false.
      do col = 1, size(model%b, 2)
         do row = 1, size(model%b, 1)
            adjacency(row, col) = potential_matrix_nonzero(constraints%b, model%b, model%b_t, row, col, t)
         end do
      end do
   end subroutine potential_b_adjacency_at

   pure subroutine kem_state_classes(adjacency, q_zero, deterministic, indirectly_stochastic)
      logical, intent(in) :: adjacency(:, :) !! Potential B adjacency matrix at one time slice.
      logical, intent(in) :: q_zero(:) !! True for states whose process-variance diagonal is fixed at zero.
      logical, intent(out) :: deterministic(:) !! True for zero-Q states with no length-m path to a directly stochastic state.
      logical, intent(out) :: indirectly_stochastic(:) !! True for zero-Q states linked through B to a positive-Q state.
      logical, allocatable :: next_reach(:, :)
      logical, allocatable :: reach(:, :)
      logical :: linked
      integer :: col
      integer :: i
      integer :: j
      integer :: power
      integer :: row

      allocate(reach(size(adjacency, 1), size(adjacency, 2)))
      allocate(next_reach(size(adjacency, 1), size(adjacency, 2)))
      reach = adjacency
      do power = 2, size(adjacency, 1)
         next_reach = .false.
         do col = 1, size(adjacency, 2)
            do row = 1, size(adjacency, 1)
               do j = 1, size(adjacency, 2)
                  if (reach(row, j) .and. adjacency(j, col)) then
                     next_reach(row, col) = .true.
                     exit
                  end if
               end do
            end do
         end do
         reach = next_reach
      end do

      deterministic = .false.
      indirectly_stochastic = .false.
      do i = 1, size(q_zero)
         if (.not. q_zero(i)) cycle
         linked = .false.
         do col = 1, size(q_zero)
            if (q_zero(col)) cycle
            if (reach(i, col)) then
               linked = .true.
               exit
            end if
         end do
         if (linked) then
            indirectly_stochastic(i) = .true.
         else
            deterministic(i) = .true.
         end if
      end do
   end subroutine kem_state_classes

   pure logical function potential_matrix_nonzero(block, base, timed, row, col, t) result(can_be_nonzero)
      type(marss_constraint_block), intent(in) :: block !! Affine constraint block used to detect potentially nonzero entries.
      real(dp), intent(in) :: base(:, :) !! Static matrix containing the current numerical value.
      real(dp), allocatable, intent(in) :: timed(:, :, :) !! Optional time-indexed matrix representation.
      integer, intent(in) :: row !! One-based matrix row.
      integer, intent(in) :: col !! One-based matrix column.
      integer, intent(in) :: t !! One-based time index used for a time-varying block.
      real(dp) :: value
      integer :: idx
      integer :: tslice

      if (allocated(timed)) then
         value = timed(row, col, t)
         tslice = t
      else
         value = base(row, col)
         tslice = 1
      end if
      idx = row + (col - 1) * size(base, 1) + (tslice - 1) * size(base)
      can_be_nonzero = abs(value) > 128.0_dp * epsilon(1.0_dp)
      if (allocated(block%fixed)) then
         if (idx >= 1 .and. idx <= size(block%fixed)) then
            can_be_nonzero = can_be_nonzero .or. abs(block%fixed(idx)) > 128.0_dp * epsilon(1.0_dp)
         end if
      end if
      if (allocated(block%design)) can_be_nonzero = can_be_nonzero .or. constraint_element_free(block, idx)
   end function potential_matrix_nonzero

   pure logical function constraint_element_free(block, idx) result(is_free)
      type(marss_constraint_block), intent(in) :: block !! Affine constraint block whose design row is inspected.
      integer, intent(in) :: idx !! One-based flattened element index in Fortran column-major order.

      is_free = .false.
      if (.not. allocated(block%design)) return
      if (idx < 1 .or. idx > size(block%design, 1)) return
      if (size(block%design, 2) == 0) return
      is_free = any(abs(block%design(idx, :)) > 128.0_dp * epsilon(1.0_dp))
   end function constraint_element_free

   pure logical function covariance_valid(a) result(ok)
      real(dp), intent(in) :: a(:, :) !! Candidate covariance matrix required to be symmetric positive semidefinite.

      ok = .false.
      if (size(a, 1) /= size(a, 2)) return
      if (.not. is_symmetric(a)) return
      if (.not. is_positive_semidefinite(a)) return
      ok = .true.
   end function covariance_valid

end module marss_model_ops
