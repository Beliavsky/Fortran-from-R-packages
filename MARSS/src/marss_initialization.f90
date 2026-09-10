! SPDX-License-Identifier: GPL-2.0-only
module marss_initialization
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_constraint_block, marss_constraints
   use marss_constraints_mod, only : marss_apply_constraints, marss_constraints_update_start_named
   use marss_parameters, only : marss_b_at, marss_u_at, marss_z_at, marss_a_at
   use marss_model_ops, only : marss_model_valid
   use marss_utils, only : identity_matrix
   use marss_covariance, only : jacobi_symmetric_eigen
   implicit none
   private
   public :: marss_inits_linear
   public :: marss_inits_named

contains

   subroutine marss_inits_linear(template, constraints, initialized, model, info)
      type(marss_model), intent(in) :: template !! Numerical model defining dimensions, data, and fixed parameter offsets.
      type(marss_constraints), intent(in) :: constraints !! Affine f+D*beta model structure whose start coordinates are initialized.
      type(marss_constraints), intent(out) :: initialized !! Constraint structure with MARSS-style default beta starting values.
      type(marss_model), intent(out) :: model !! Numerical model obtained by applying the initialized beta coordinates.
      integer, intent(out) :: info !! Zero on success; positive values identify projection, x0-solve, or model-validation failures.
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: target(:)
      integer :: p
      integer :: slice_size

      initialized = constraints
      info = 0

      call repeated_diagonal_target(template, 'b', 1.0_dp, target, info, slice_size)
      if (info /= 0) then
         info = 10 + info
         return
      end if
      call project_start(initialized%b, target, info, slice_size)
      if (info /= 0) then
         info = 20 + info
         return
      end if

      call fill_start(initialized%u, 0.0_dp)

      call repeated_diagonal_target(template, 'q', 0.05_dp, target, info, slice_size)
      if (info /= 0) then
         info = 30 + info
         return
      end if
      call project_start(initialized%q, target, info, slice_size)
      if (info /= 0) then
         info = 40 + info
         return
      end if

      call fill_start(initialized%z, 1.0_dp)
      call fill_start(initialized%a, 0.0_dp)

      call repeated_diagonal_target(template, 'r', 0.05_dp, target, info, slice_size)
      if (info /= 0) then
         info = 50 + info
         return
      end if
      call project_start(initialized%r, target, info, slice_size)
      if (info /= 0) then
         info = 60 + info
         return
      end if

      call repeated_diagonal_target(template, 'v', 5.0_dp, target, info, slice_size)
      if (info /= 0) then
         info = 70 + info
         return
      end if
      call project_start(initialized%v0, target, info, slice_size)
      if (info /= 0) then
         info = 80 + info
         return
      end if

      call fill_start(initialized%c, 0.0_dp)
      call fill_start(initialized%d, 0.0_dp)
      call fill_start(initialized%g, 0.0_dp)
      call fill_start(initialized%h, 0.0_dp)
      call fill_start(initialized%l, 0.0_dp)
      call fill_start(initialized%x0, 0.0_dp)

      call constraints_start_vector(initialized, theta)
      call marss_apply_constraints(template, initialized, theta, model, info)
      if (info /= 0) then
         info = 90 + info
         return
      end if

      p = block_parameter_count(initialized%x0)
      if (p > 0) then
         call initialize_x0(model, initialized%x0, info)
         if (info /= 0) then
            info = 100 + info
            return
         end if
         call constraints_start_vector(initialized, theta)
         call marss_apply_constraints(template, initialized, theta, model, info)
         if (info /= 0) then
            info = 110 + info
            return
         end if
      end if

      if (.not. marss_model_valid(model)) then
         info = 120
         return
      end if
   end subroutine marss_inits_linear

   subroutine marss_inits_named(template, constraints, parameter_names, values, initialized, model, info)
      type(marss_model), intent(in) :: template !! Numerical model defining data, dimensions, and fixed parameter offsets.
      type(marss_constraints), intent(in) :: constraints !! Affine model structure to initialize with MARSS defaults.
      character(len=*), intent(in) :: parameter_names(:) !! Raw or block-prefixed free-coordinate names to override.
      real(dp), intent(in) :: values(:) !! User starting values for the named subset of free beta coordinates.
      type(marss_constraints), intent(out) :: initialized !! Default-initialized constraints with named overrides applied.
      type(marss_model), intent(out) :: model !! Numerical model obtained from the overridden initialization coordinates.
      integer, intent(out) :: info !! Zero on success; nonzero for default initialization, name matching, or model validity failure.
      type(marss_constraints) :: defaults
      type(marss_model) :: default_model
      real(dp), allocatable :: theta(:)
      integer :: local_info

      call marss_inits_linear(template, constraints, defaults, default_model, local_info)
      if (local_info /= 0) then
         info = local_info
         return
      end if
      call marss_constraints_update_start_named(defaults, parameter_names, values, initialized, local_info)
      if (local_info /= 0) then
         info = 200 + local_info
         return
      end if
      call constraints_start_vector(initialized, theta)
      call marss_apply_constraints(template, initialized, theta, model, local_info)
      if (local_info /= 0) then
         info = 210 + local_info
         return
      end if
      if (.not. marss_model_valid(model)) then
         info = 220
         return
      end if
      info = 0
   end subroutine marss_inits_named

   pure subroutine fill_start(block, value)
      type(marss_constraint_block), intent(inout) :: block !! Constraint block whose beta start coordinates are filled uniformly.
      real(dp), intent(in) :: value !! Scalar value assigned to every free beta coordinate in the block.

      if (.not. allocated(block%start)) return
      block%start = value
   end subroutine fill_start

   pure integer function block_parameter_count(block) result(nparam)
      type(marss_constraint_block), intent(in) :: block !! Constraint block whose number of beta coordinates is requested.

      nparam = 0
      if (allocated(block%design)) nparam = size(block%design, 2)
   end function block_parameter_count

   pure subroutine repeated_diagonal_target(model, which, diagonal, target, info, slice_size)
      type(marss_model), intent(in) :: model !! Model defining raw parameter dimensions and optional time indexing.
      character(len=1), intent(in) :: which !! Block selector: b, q, r, or v for B, Q, R, or V0.
      real(dp), intent(in) :: diagonal !! Scalar diagonal value used by MARSS default initialization.
      real(dp), allocatable, intent(out) :: target(:) !! Flattened target matrix or repeated time-indexed matrices.
      integer, intent(out) :: info !! Zero on success, one for an unknown selector.
      integer, intent(out) :: slice_size !! Number of flattened entries in one time slice of the selected matrix.
      real(dp), allocatable :: matrix(:, :)
      integer :: n
      integer :: nslice
      integer :: t

      info = 0
      select case (which)
      case ('b')
         n = size(model%b, 1)
         nslice = 1
         if (allocated(model%b_t)) nslice = size(model%b_t, 3)
      case ('q')
         if (allocated(model%q_noise)) then
            n = size(model%q_noise, 1)
            nslice = 1
            if (allocated(model%q_noise_t)) nslice = size(model%q_noise_t, 3)
         else
            n = size(model%q, 1)
            nslice = 1
            if (allocated(model%q_t)) nslice = size(model%q_t, 3)
         end if
      case ('r')
         if (allocated(model%r_noise)) then
            n = size(model%r_noise, 1)
            nslice = 1
            if (allocated(model%r_noise_t)) nslice = size(model%r_noise_t, 3)
         else
            n = size(model%r, 1)
            nslice = 1
            if (allocated(model%r_t)) nslice = size(model%r_t, 3)
         end if
      case ('v')
         if (allocated(model%v0_noise)) then
            n = size(model%v0_noise, 1)
         else
            n = size(model%v0, 1)
         end if
         nslice = 1
      case default
         info = 1
         slice_size = 0
         allocate(target(0))
         return
      end select
      slice_size = n * n
      matrix = diagonal * identity_matrix(n)
      allocate(target(n * n * nslice))
      do t = 1, nslice
         target((t - 1) * n * n + 1:t * n * n) = reshape(matrix, [n * n])
      end do
   end subroutine repeated_diagonal_target

   pure subroutine project_start(block, target, info, slice_size)
      type(marss_constraint_block), intent(inout) :: block !! Affine block whose default target is projected into beta coordinates.
      real(dp), intent(in) :: target(:) !! Desired numerical values before respecting fixed/equality constraints.
      integer, intent(out) :: info !! Zero on success; nonzero for incompatible dimensions or an eigensolver failure.
      integer, intent(in), optional :: slice_size !! Raw entries per time slice; enables upstream time-averaged scalar projection.
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: averaged_design(:, :)
      real(dp), allocatable :: averaged_fixed(:)
      real(dp), allocatable :: desired(:)
      real(dp), allocatable :: gram(:, :)
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: eigenvectors(:, :)
      real(dp), allocatable :: work(:)
      real(dp) :: tolerance
      integer :: count_nonzero
      integer :: eig_info
      integer :: i
      integer :: j
      integer :: nslice
      integer :: p
      integer :: t

      info = 0
      if (.not. allocated(block%design)) return
      if (.not. allocated(block%fixed) .or. .not. allocated(block%start)) then
         info = 1
         return
      end if
      if (size(target) /= size(block%fixed)) then
         info = 2
         return
      end if
      p = size(block%design, 2)
      if (p == 0) return
      if (present(slice_size) .and. slice_size > 0 .and. size(target) > slice_size) then
         if (mod(size(target), slice_size) /= 0) then
            info = 4
            return
         end if
         nslice = size(target) / slice_size
         allocate(averaged_design(slice_size, p), averaged_fixed(slice_size), desired(slice_size))
         averaged_design = 0.0_dp
         averaged_fixed = 0.0_dp
         desired = target(1:slice_size)
         do i = 1, slice_size
            if (all(abs(block%design(i:size(block%fixed):slice_size, :)) < tiny(1.0_dp))) then
               desired(i) = block%fixed(i)
            end if
            count_nonzero = 0
            do t = 1, nslice
               j = (t - 1) * slice_size + i
               if (abs(block%fixed(j)) > tiny(1.0_dp)) then
                  averaged_fixed(i) = averaged_fixed(i) + block%fixed(j)
                  count_nonzero = count_nonzero + 1
               end if
            end do
            if (count_nonzero > 0) averaged_fixed(i) = averaged_fixed(i) / real(count_nonzero, dp)
            do j = 1, p
               count_nonzero = 0
               do t = 1, nslice
                  if (abs(block%design((t - 1) * slice_size + i, j)) > tiny(1.0_dp)) then
                     averaged_design(i, j) = averaged_design(i, j) + &
                        block%design((t - 1) * slice_size + i, j)
                     count_nonzero = count_nonzero + 1
                  end if
               end do
               if (count_nonzero > 0) averaged_design(i, j) = averaged_design(i, j) / real(count_nonzero, dp)
            end do
         end do
         rhs = matmul(transpose(averaged_design), desired - averaged_fixed)
         gram = matmul(transpose(averaged_design), averaged_design)
      else
         rhs = matmul(transpose(block%design), target - block%fixed)
         gram = matmul(transpose(block%design), block%design)
      end if
      call jacobi_symmetric_eigen(gram, eigenvalues, eigenvectors, eig_info)
      if (eig_info /= 0) then
         info = 3
         return
      end if
      work = matmul(transpose(eigenvectors), rhs)
      tolerance = 128.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(eigenvalues)))
      do i = 1, p
         if (abs(eigenvalues(i)) > tolerance) then
            work(i) = work(i) / eigenvalues(i)
         else
            work(i) = 0.0_dp
         end if
      end do
      block%start = matmul(eigenvectors, work)
   end subroutine project_start

   subroutine initialize_x0(model, block, info)
      type(marss_model), intent(in) :: model !! Model with all non-x0 default starts applied before solving the initial state.
      type(marss_constraint_block), intent(inout) :: block !! x0 affine block whose beta start is solved from the first observation.
      integer, intent(out) :: info !! Zero on success; nonzero when x0 is underconstrained or the solve fails.
      real(dp), allocatable :: design_observation(:, :)
      real(dp), allocatable :: gram(:, :)
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: eigenvectors(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: work(:)
      real(dp) :: b1(size(model%b, 1), size(model%b, 2))
      real(dp) :: z1(size(model%z, 1), size(model%z, 2))
      real(dp) :: u1(size(model%u))
      real(dp) :: a1(size(model%a))
      real(dp) :: y1(size(model%y, 1))
      real(dp) :: tolerance
      integer :: eig_info
      integer :: i
      integer :: p
      integer :: rank

      info = 0
      p = block_parameter_count(block)
      if (p == 0) return
      if (.not. allocated(block%fixed) .or. .not. allocated(block%design)) then
         info = 1
         return
      end if
      y1 = model%y(:, 1)
      do i = 1, size(y1)
         if (ieee_is_nan(y1(i))) y1(i) = 0.0_dp
      end do
      z1 = marss_z_at(model, 1)
      a1 = marss_a_at(model, 1)
      if (model%tinitx == 0) then
         b1 = marss_b_at(model, 1)
         u1 = marss_u_at(model, 1)
         design_observation = matmul(matmul(z1, b1), block%design)
         rhs = y1 - matmul(z1, matmul(b1, block%fixed) + u1) - a1
      else
         design_observation = matmul(z1, block%design)
         rhs = y1 - matmul(z1, block%fixed) - a1
      end if
      gram = matmul(transpose(design_observation), design_observation)
      call jacobi_symmetric_eigen(gram, eigenvalues, eigenvectors, eig_info)
      if (eig_info /= 0) then
         info = 2
         return
      end if
      tolerance = 128.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(eigenvalues)))
      rank = count(abs(eigenvalues) > tolerance)
      if (rank < p) then
         info = 3
         return
      end if
      work = matmul(transpose(eigenvectors), matmul(transpose(design_observation), rhs))
      do i = 1, p
         if (abs(eigenvalues(i)) > tolerance) then
            work(i) = work(i) / eigenvalues(i)
         else
            work(i) = 0.0_dp
         end if
      end do
      block%start = matmul(eigenvectors, work)
   end subroutine initialize_x0

   pure subroutine constraints_start_vector(constraints, theta)
      type(marss_constraints), intent(in) :: constraints !! Affine blocks whose current beta starts are concatenated.
      real(dp), allocatable, intent(out) :: theta(:) !! Beta coordinates in B,U,Q,Z,A,R,x0,V0,C,D,G,H,L order.
      integer :: k
      integer :: n

      n = block_parameter_count(constraints%b) + block_parameter_count(constraints%u) + &
         block_parameter_count(constraints%q) + block_parameter_count(constraints%z) + &
         block_parameter_count(constraints%a) + block_parameter_count(constraints%r) + &
         block_parameter_count(constraints%x0) + block_parameter_count(constraints%v0) + &
         block_parameter_count(constraints%c) + block_parameter_count(constraints%d) + &
         block_parameter_count(constraints%g) + block_parameter_count(constraints%h) + &
         block_parameter_count(constraints%l)
      allocate(theta(n))
      k = 0
      call append_start(constraints%b, theta, k)
      call append_start(constraints%u, theta, k)
      call append_start(constraints%q, theta, k)
      call append_start(constraints%z, theta, k)
      call append_start(constraints%a, theta, k)
      call append_start(constraints%r, theta, k)
      call append_start(constraints%x0, theta, k)
      call append_start(constraints%v0, theta, k)
      call append_start(constraints%c, theta, k)
      call append_start(constraints%d, theta, k)
      call append_start(constraints%g, theta, k)
      call append_start(constraints%h, theta, k)
      call append_start(constraints%l, theta, k)
   end subroutine constraints_start_vector

   pure subroutine append_start(block, theta, k)
      type(marss_constraint_block), intent(in) :: block !! Constraint block whose current beta start is appended.
      real(dp), intent(inout) :: theta(:) !! Destination concatenated beta vector.
      integer, intent(inout) :: k !! Number of coordinates already written, updated on return.
      integer :: p

      if (.not. allocated(block%start)) return
      p = size(block%start)
      if (p > 0) theta(k + 1:k + p) = block%start
      k = k + p
   end subroutine append_start

end module marss_initialization
