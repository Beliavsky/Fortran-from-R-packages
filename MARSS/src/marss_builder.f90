! SPDX-License-Identifier: GPL-2.0-only
module marss_builder
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_model_spec, marss_constraint_block, marss_constraints
   use marss_utils, only : identity_matrix
   use marss_constraints_mod, only : marss_apply_constraints, marss_constraint_start_vector
   use marss_model_ops, only : marss_model_valid
   implicit none
   private
   public :: marss_build
   public :: canonical_spec
   public :: make_fixed_block
   public :: make_fixed_vector
   public :: make_vector_block
   public :: make_matrix_block

contains

   subroutine marss_build(y, spec, model, constraints, info, state_covariates, obs_covariates)
      real(dp), intent(in) :: y(:, :) !! Observation matrix, variables by time; NaNs denote missing observations.
      type(marss_model_spec), intent(in) :: spec !! High-level MARSS shortcut specification and state dimension.
      type(marss_model), intent(out) :: model !! Numerical model initialized consistently with the requested shortcuts.
      type(marss_constraints), intent(out) :: constraints !! Affine fixed/free/equality representation of the requested shortcuts.
      integer, intent(out) :: info !! Zero on success; positive values identify an invalid or unsupported shortcut combination.
      real(dp), intent(in), optional :: state_covariates(:, :) !! Optional c(t) rows; one column is recycled over time.
      real(dp), intent(in), optional :: obs_covariates(:, :) !! Optional d(t) rows; one column is recycled over time.
      type(marss_model) :: template
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: target(:, :)
      real(dp), allocatable :: target_vec(:)
      integer :: i
      integer :: m
      integer :: n
      integer :: tt
      character(len=32) :: zspec
      character(len=32) :: form

      info = 0
      n = size(y, 1)
      tt = size(y, 2)
      if (n < 1 .or. tt < 1) then
         info = 1
         return
      end if
      form = canonical_spec(spec%form)
      if (form /= "marss" .and. form /= "marxss") then
         info = 5
         return
      end if
      zspec = canonical_spec(spec%z)
      m = spec%nstate
      if (m <= 0) then
         if (zspec == "onestate") then
            m = 1
         else
            m = n
         end if
      end if
      if (m < 1) then
         info = 2
         return
      end if
      if (spec%tinitx /= 0 .and. spec%tinitx /= 1) then
         info = 3
         return
      end if
      if (spec%diffuse .and. spec%tinitx /= 1) then
         info = 4
         return
      end if

      allocate(template%y(n, tt), template%b(m, m), template%u(m), template%q(m, m))
      allocate(template%z(n, m), template%a(n), template%r(n, n))
      allocate(template%x0(m), template%v0(m, m))
      allocate(template%g(m, m), template%h(n, n), template%l(m, m))
      allocate(template%q_noise(m, m), template%r_noise(n, n), template%v0_noise(m, m))
      template%y = y
      template%b = identity_matrix(m)
      template%u = 0.0_dp
      template%q = 0.05_dp * identity_matrix(m)
      template%z = initial_z(n, m, zspec)
      template%a = 0.0_dp
      template%r = 0.05_dp * identity_matrix(n)
      template%x0 = initial_x0(y, m)
      template%v0 = 5.0_dp * identity_matrix(m)
      template%q_noise = template%q
      template%r_noise = template%r
      template%v0_noise = template%v0
      template%tinitx = spec%tinitx
      template%diffuse = spec%diffuse

      call fixed_loading(canonical_spec(spec%g), m, template%g, constraints%g, info)
      if (info /= 0) then
         info = 100 + info
         return
      end if
      call fixed_loading(canonical_spec(spec%h), n, template%h, constraints%h, info)
      if (info /= 0) then
         info = 110 + info
         return
      end if
      call fixed_loading(canonical_spec(spec%l), m, template%l, constraints%l, info)
      if (info /= 0) then
         info = 120 + info
         return
      end if

      target = identity_matrix(m)
      call make_matrix_block(canonical_spec(spec%b), target, .false., constraints%b, info)
      if (info /= 0) then
         info = 200 + info
         return
      end if
      target_vec = spread(0.0_dp, 1, m)
      call make_vector_block(canonical_spec(spec%u), target_vec, constraints%u, info)
      if (info /= 0) then
         info = 210 + info
         return
      end if
      target = 0.05_dp * identity_matrix(m)
      call make_matrix_block(canonical_spec(spec%q), target, .true., constraints%q, info)
      if (info /= 0) then
         info = 220 + info
         return
      end if

      if (zspec == "onestate") then
         if (m /= 1) then
            info = 230
            return
         end if
         target = reshape([(1.0_dp, i = 1, n)], [n, 1])
         call make_fixed_block(target, constraints%z)
      else
         target = initial_z(n, m, zspec)
         call make_matrix_block(zspec, target, .false., constraints%z, info)
         if (info /= 0) then
            info = 231 + info
            return
         end if
      end if

      target_vec = spread(0.0_dp, 1, n)
      if (canonical_spec(spec%a) == "scaling") then
         call make_scaling_a(template%z, constraints%z, target_vec, constraints%a, info)
      else
         call make_vector_block(canonical_spec(spec%a), target_vec, constraints%a, info)
      end if
      if (info /= 0) then
         info = 240 + info
         return
      end if

      target = 0.05_dp * identity_matrix(n)
      call make_matrix_block(canonical_spec(spec%r), target, .true., constraints%r, info)
      if (info /= 0) then
         info = 250 + info
         return
      end if

      target_vec = initial_x0(y, m)
      call make_vector_block(canonical_spec(spec%x0), target_vec, constraints%x0, info)
      if (info /= 0) then
         info = 260 + info
         return
      end if
      target = 5.0_dp * identity_matrix(m)
      call make_matrix_block(canonical_spec(spec%v0), target, .true., constraints%v0, info)
      if (info /= 0) then
         info = 270 + info
         return
      end if

      if (present(state_covariates)) then
         call install_covariates(state_covariates, tt, template%state_covariates, info)
         if (info /= 0) then
            info = 280 + info
            return
         end if
         allocate(template%c_coef(m, size(template%state_covariates, 1)))
         template%c_coef = 0.0_dp
         target = template%c_coef
         call make_matrix_block(canonical_spec(spec%c), target, .false., constraints%c, info)
         if (info /= 0) then
            info = 281 + info
            return
         end if
      else if (canonical_spec(spec%c) /= "zero") then
         info = 282
         return
      end if

      if (present(obs_covariates)) then
         call install_covariates(obs_covariates, tt, template%obs_covariates, info)
         if (info /= 0) then
            info = 290 + info
            return
         end if
         allocate(template%d_coef(n, size(template%obs_covariates, 1)))
         template%d_coef = 0.0_dp
         target = template%d_coef
         call make_matrix_block(canonical_spec(spec%d), target, .false., constraints%d, info)
         if (info /= 0) then
            info = 291 + info
            return
         end if
      else if (canonical_spec(spec%d) /= "zero") then
         info = 292
         return
      end if

      call marss_constraint_start_vector(constraints, theta)
      call marss_apply_constraints(template, constraints, theta, model, info)
      if (info /= 0) then
         info = 300 + info
         return
      end if

      model%q = matmul(matmul(model%g, model%q_noise), transpose(model%g))
      model%r = matmul(matmul(model%h, model%r_noise), transpose(model%h))
      model%v0 = matmul(matmul(model%l, model%v0_noise), transpose(model%l))
      if (.not. marss_model_valid(model)) then
         info = 400
         return
      end if

      do i = 1, m
         if (ieee_is_nan(model%x0(i))) model%x0(i) = 0.0_dp
      end do
   end subroutine marss_build

   pure function canonical_spec(text) result(value)
      character(len=*), intent(in) :: text !! Model shortcut text to normalize to lowercase with collapsed surrounding blanks.
      character(len=32) :: value
      integer :: i
      integer :: code

      value = ""
      value = adjustl(trim(text))
      do i = 1, len_trim(value)
         code = iachar(value(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) value(i:i) = achar(code + 32)
      end do
   end function canonical_spec

   pure function initial_z(n, m, spec) result(z)
      integer, intent(in) :: n !! Number of observation series.
      integer, intent(in) :: m !! Number of latent states.
      character(len=*), intent(in) :: spec !! Canonical Z shortcut used to select a stable numerical starting matrix.
      real(dp) :: z(n, m)
      integer :: i

      z = 0.0_dp
      if (spec == "onestate") then
         z = 1.0_dp
      else if (spec == "unconstrained" .or. spec == "unequal" .or. spec == "equal") then
         z = 1.0_dp
      else
         do i = 1, min(n, m)
            z(i, i) = 1.0_dp
         end do
      end if
   end function initial_z

   pure function initial_x0(y, m) result(x0)
      real(dp), intent(in) :: y(:, :) !! Observation matrix used to seed initial state values from the first available observations.
      integer, intent(in) :: m !! Number of state components to initialize.
      real(dp) :: x0(m)
      integer :: i
      integer :: t

      x0 = 0.0_dp
      do i = 1, min(m, size(y, 1))
         do t = 1, size(y, 2)
            if (.not. ieee_is_nan(y(i, t))) then
               x0(i) = y(i, t)
               exit
            end if
         end do
      end do
   end function initial_x0

   pure subroutine install_covariates(input, tt, output, info)
      real(dp), intent(in) :: input(:, :) !! Covariate matrix with one or TT columns.
      integer, intent(in) :: tt !! Number of observation times required by the numerical model.
      real(dp), allocatable, intent(out) :: output(:, :) !! Covariates expanded to one column per observation time.
      integer, intent(out) :: info !! Zero on success, one when the time dimension is neither one nor TT.
      integer :: t

      info = 0
      if (size(input, 1) < 1) then
         info = 1
         return
      end if
      if (size(input, 2) == tt) then
         output = input
      else if (size(input, 2) == 1) then
         allocate(output(size(input, 1), tt))
         do t = 1, tt
            output(:, t) = input(:, 1)
         end do
      else
         info = 1
      end if
   end subroutine install_covariates

   pure subroutine fixed_loading(spec, n, loading, block, info)
      character(len=*), intent(in) :: spec !! Canonical G, H, or L shortcut; only identity and zero are valid in form=marss.
      integer, intent(in) :: n !! Square loading dimension.
      real(dp), intent(out) :: loading(:, :) !! Fixed loading matrix represented by the shortcut.
      type(marss_constraint_block), intent(out) :: block !! Zero-parameter affine block recording the fixed loading.
      integer, intent(out) :: info !! Zero on success, one for an unsupported loading shortcut.

      info = 0
      select case (spec)
      case ("identity")
         loading = identity_matrix(n)
      case ("zero")
         loading = 0.0_dp
      case default
         info = 1
         return
      end select
      call make_fixed_block(loading, block)
   end subroutine fixed_loading

   pure subroutine make_fixed_block(values, block)
      real(dp), intent(in) :: values(:, :) !! Fixed matrix values in ordinary row/column form.
      type(marss_constraint_block), intent(out) :: block !! Affine block with zero free coordinates.

      allocate(block%fixed(size(values)), block%design(size(values), 0), block%start(0))
      block%fixed = reshape(values, [size(values)])
   end subroutine make_fixed_block

   pure subroutine make_fixed_vector(values, block)
      real(dp), intent(in) :: values(:) !! Fixed vector values.
      type(marss_constraint_block), intent(out) :: block !! Affine block with zero free coordinates.

      allocate(block%fixed(size(values)), block%design(size(values), 0), block%start(0))
      block%fixed = values
   end subroutine make_fixed_vector

   pure subroutine make_vector_block(spec, target, block, info)
      character(len=*), intent(in) :: spec !! Canonical vector shortcut: unconstrained/unequal, equal, zero, or diagonal aliases.
      real(dp), intent(in) :: target(:) !! Desired numerical starting vector used to initialize beta coordinates.
      type(marss_constraint_block), intent(out) :: block !! Affine representation of the vector shortcut.
      integer, intent(out) :: info !! Zero on success, one for an unsupported vector shortcut.
      integer :: i
      integer :: n

      info = 0
      n = size(target)
      select case (spec)
      case ("zero")
         call make_fixed_vector(0.0_dp * target, block)
      case ("unconstrained", "unequal", "diagonal and unequal")
         allocate(block%fixed(n), block%design(n, n), block%start(n))
         block%fixed = 0.0_dp
         block%design = 0.0_dp
         do i = 1, n
            block%design(i, i) = 1.0_dp
         end do
         block%start = target
      case ("equal", "diagonal and equal")
         allocate(block%fixed(n), block%design(n, 1), block%start(1))
         block%fixed = 0.0_dp
         block%design = 1.0_dp
         if (n > 0) then
            block%start(1) = sum(target) / real(n, dp)
         else
            block%start(1) = 0.0_dp
         end if
      case default
         info = 1
      end select
   end subroutine make_vector_block

   pure subroutine make_matrix_block(spec, target, symmetric, block, info)
      character(len=*), intent(in) :: spec !! Canonical matrix shortcut from the MARSS form=marss model vocabulary.
      real(dp), intent(in) :: target(:, :) !! Desired numerical starting matrix used to initialize free coordinates.
      logical, intent(in) :: symmetric !! True when free off-diagonal entries must share covariance symmetry parameters.
      type(marss_constraint_block), intent(out) :: block !! Affine representation of the requested shortcut.
      integer, intent(out) :: info !! Zero on success; positive values indicate a shape or shortcut error.
      integer :: col
      integer :: i
      integer :: idx
      integer :: j
      integer :: ncol
      integer :: nelem
      integer :: nrow
      integer :: p
      integer :: row

      info = 0
      nrow = size(target, 1)
      ncol = size(target, 2)
      nelem = size(target)
      select case (spec)
      case ("identity")
         if (nrow /= ncol) then
            info = 1
            return
         end if
         call make_fixed_block(identity_matrix(nrow), block)
      case ("zero")
         call make_fixed_block(0.0_dp * target, block)
      case ("unconstrained", "unequal")
         if (symmetric) then
            if (nrow /= ncol) then
               info = 2
               return
            end if
            p = nrow * (nrow + 1) / 2
            allocate(block%fixed(nelem), block%design(nelem, p), block%start(p))
            block%fixed = 0.0_dp
            block%design = 0.0_dp
            p = 0
            do col = 1, ncol
               do row = col, nrow
                  p = p + 1
                  idx = row + (col - 1) * nrow
                  block%design(idx, p) = 1.0_dp
                  idx = col + (row - 1) * nrow
                  block%design(idx, p) = 1.0_dp
                  block%start(p) = target(row, col)
               end do
            end do
         else
            allocate(block%fixed(nelem), block%design(nelem, nelem), block%start(nelem))
            block%fixed = 0.0_dp
            block%design = 0.0_dp
            do i = 1, nelem
               block%design(i, i) = 1.0_dp
            end do
            block%start = reshape(target, [nelem])
         end if
      case ("diagonal and unequal")
         if (nrow /= ncol) then
            info = 3
            return
         end if
         p = nrow
         allocate(block%fixed(nelem), block%design(nelem, p), block%start(p))
         block%fixed = 0.0_dp
         block%design = 0.0_dp
         do i = 1, nrow
            idx = i + (i - 1) * nrow
            block%design(idx, i) = 1.0_dp
            block%start(i) = target(i, i)
         end do
      case ("diagonal and equal")
         if (nrow /= ncol) then
            info = 4
            return
         end if
         allocate(block%fixed(nelem), block%design(nelem, 1), block%start(1))
         block%fixed = 0.0_dp
         block%design = 0.0_dp
         do i = 1, nrow
            idx = i + (i - 1) * nrow
            block%design(idx, 1) = 1.0_dp
         end do
         block%start(1) = sum([(target(i, i), i = 1, nrow)]) / real(nrow, dp)
      case ("equalvarcov")
         if (nrow /= ncol) then
            info = 5
            return
         end if
         p = merge(1, 2, nrow == 1)
         allocate(block%fixed(nelem), block%design(nelem, p), block%start(p))
         block%fixed = 0.0_dp
         block%design = 0.0_dp
         do col = 1, ncol
            do row = 1, nrow
               idx = row + (col - 1) * nrow
               if (row == col) then
                  block%design(idx, 1) = 1.0_dp
               else if (p == 2) then
                  block%design(idx, 2) = 1.0_dp
               end if
            end do
         end do
         block%start(1) = sum([(target(i, i), i = 1, nrow)]) / real(nrow, dp)
         if (p == 2) then
            block%start(2) = 0.0_dp
            p = 0
            do j = 1, ncol
               do i = 1, nrow
                  if (i == j) cycle
                  block%start(2) = block%start(2) + target(i, j)
                  p = p + 1
               end do
            end do
            if (p > 0) block%start(2) = block%start(2) / real(p, dp)
         end if
      case ("equal")
         allocate(block%fixed(nelem), block%design(nelem, 1), block%start(1))
         block%fixed = 0.0_dp
         block%design = 1.0_dp
         block%start(1) = sum(target) / real(nelem, dp)
      case default
         info = 6
      end select
   end subroutine make_matrix_block

   pure subroutine make_scaling_a(z, zblock, target, block, info)
      real(dp), intent(in) :: z(:, :) !! Numerical time-constant observation design matrix.
      type(marss_constraint_block), intent(in) :: zblock !! Z constraint block, which must be completely fixed for scaling A.
      real(dp), intent(in) :: target(:) !! Desired starting observation intercepts, normally zero.
      type(marss_constraint_block), intent(out) :: block !! A block with one fixed zero reference per state-loading column.
      integer, intent(out) :: info !! Zero on success; positive when Z is not a fixed 0/1 design matrix.
      logical, allocatable :: reference(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: nfree
      real(dp) :: tolerance

      info = 0
      if (.not. allocated(zblock%design)) then
         info = 1
         return
      end if
      if (size(zblock%design, 2) /= 0) then
         info = 2
         return
      end if
      tolerance = 128.0_dp * epsilon(1.0_dp)
      do i = 1, size(z, 1)
         if (any(abs(z(i, :)) > tolerance .and. abs(z(i, :) - 1.0_dp) > tolerance)) then
            info = 3
            return
         end if
         if (abs(sum(z(i, :)) - 1.0_dp) > tolerance) then
            info = 4
            return
         end if
      end do
      allocate(reference(size(z, 1)))
      reference = .false.
      do j = 1, size(z, 2)
         do i = 1, size(z, 1)
            if (abs(z(i, j)) > tolerance) then
               reference(i) = .true.
               exit
            end if
         end do
      end do
      nfree = count(.not. reference)
      allocate(block%fixed(size(target)), block%design(size(target), nfree), block%start(nfree))
      block%fixed = 0.0_dp
      block%design = 0.0_dp
      block%start = 0.0_dp
      k = 0
      do i = 1, size(target)
         if (reference(i)) cycle
         k = k + 1
         block%design(i, k) = 1.0_dp
         block%start(k) = target(i)
      end do
   end subroutine make_scaling_a

end module marss_builder
