! SPDX-License-Identifier: GPL-2.0-only
module marss_constraints_mod
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_fit_result, marss_kf_result
   use marss_types, only : marss_constraint_block, marss_constraints, marss_parameter_name_len
   use marss_kalman, only : marss_kf
   use marss_model_ops, only : marss_model_valid
   use marss_utils, only : identity_matrix, normal_quantile
   use r_linalg, only : inverse_matrix
   implicit none
   private
   public :: marss_constraints_valid
   public :: marss_constraints_parameter_count
   public :: marss_constraint_from_labels
   public :: marss_constraint_from_affine
   public :: marss_constraint_from_entries
   public :: marss_free_parameter_names
   public :: marss_vectorized_parameter_names
   public :: marss_reorder_free_parameters
   public :: marss_constraints_set_start_named
   public :: marss_constraints_update_start_named
   public :: marss_constraint_start_vector
   public :: marss_constraints_set_start
   public :: marss_vectorize_free
   public :: marss_unvectorize_free
   public :: marss_apply_constraints
   public :: marss_optim_linear
   public :: marss_hessian_linear
   public :: marss_fisher_i_linear
   public :: marss_param_cis_linear

contains

   pure logical function marss_constraints_valid(template, constraints) result(ok)
      type(marss_model), intent(in) :: template !! Model defining the active static or time-indexed shape of each constrained block.
      type(marss_constraints), intent(in) :: constraints !! Affine f+D*beta specifications; unspecified blocks remain fixed.
      type(marss_model) :: model
      real(dp), allocatable :: theta(:)
      integer :: info

      ok = .false.
      if (.not. marss_model_valid(template)) return
      if (.not. block_valid(constraints%b, active_b_size(template))) return
      if (.not. block_valid(constraints%u, active_u_size(template))) return
      if (.not. block_valid(constraints%q, active_q_size(template))) return
      if (.not. block_valid(constraints%z, active_z_size(template))) return
      if (.not. block_valid(constraints%a, active_a_size(template))) return
      if (.not. block_valid(constraints%r, active_r_size(template))) return
      if (.not. block_valid(constraints%x0, size(template%x0))) return
      if (.not. block_valid(constraints%v0, active_v0_size(template))) return
      if (.not. block_valid(constraints%c, active_c_size(template))) return
      if (.not. block_valid(constraints%d, active_d_size(template))) return
      if (.not. block_valid(constraints%g, active_g_size(template))) return
      if (.not. block_valid(constraints%h, active_h_size(template))) return
      if (.not. block_valid(constraints%l, active_l_size(template))) return
      call marss_constraint_start_vector(constraints, theta)
      call marss_apply_constraints(template, constraints, theta, model, info)
      ok = info == 0
   end function marss_constraints_valid

   pure integer function marss_constraints_parameter_count(constraints) result(nparam)
      type(marss_constraints), intent(in) :: constraints !! Affine constraint blocks whose free-coordinate count is requested.

      nparam = block_parameter_count(constraints%b) + block_parameter_count(constraints%u) + &
         block_parameter_count(constraints%q) + block_parameter_count(constraints%z) + &
         block_parameter_count(constraints%a) + block_parameter_count(constraints%r) + &
         block_parameter_count(constraints%x0) + block_parameter_count(constraints%v0) + &
         block_parameter_count(constraints%c) + block_parameter_count(constraints%d) + &
         block_parameter_count(constraints%g) + block_parameter_count(constraints%h) + &
         block_parameter_count(constraints%l)
   end function marss_constraints_parameter_count

   pure subroutine marss_constraint_from_affine(fixed, design, start, block, info, parameter_names)
      real(dp), intent(in) :: fixed(:) !! Affine intercept f for the flattened model block.
      real(dp), intent(in) :: design(:, :) !! Affine design D mapping free coordinates into flattened model values.
      real(dp), intent(in) :: start(:) !! Starting values for the free beta coordinates.
      type(marss_constraint_block), intent(out) :: block !! Validated affine block storing f, D, beta starts, and optional names.
      integer, intent(out) :: info !! Zero on success; nonzero for incompatible dimensions or invalid parameter names.
      character(len=*), intent(in), optional :: parameter_names(:) !! Optional unique names for the free beta coordinates.
      integer :: i

      info = 0
      if (size(design, 1) /= size(fixed) .or. size(design, 2) /= size(start)) then
         info = 1
         return
      end if
      if (present(parameter_names)) then
         if (size(parameter_names) /= size(start)) then
            info = 2
            return
         end if
         do i = 1, size(parameter_names)
            if (len_trim(parameter_names(i)) == 0) then
               info = 3
               return
            end if
            if (i > 1) then
               if (any(parameter_names(1:i - 1) == parameter_names(i))) then
                  info = 4
                  return
               end if
            end if
         end do
      end if
      allocate(block%fixed(size(fixed)), block%design(size(design, 1), size(design, 2)))
      allocate(block%start(size(start)))
      block%fixed = fixed
      block%design = design
      block%start = start
      if (present(parameter_names)) then
         allocate(block%free_names(size(start)))
         block%free_names = parameter_names
      end if
   end subroutine marss_constraint_from_affine

   subroutine marss_constraint_from_entries(numeric_values, expressions, is_expression, block, info, parameter_starts)
      real(dp), intent(in) :: numeric_values(:) !! Fixed numeric entries, used where is_expression is false.
      character(len=*), intent(in) :: expressions(:) !! Character labels or affine expressions for entries marked true.
      logical, intent(in) :: is_expression(:) !! True for character/list-matrix entries; false for fixed numeric entries.
      type(marss_constraint_block), intent(out) :: block !! Affine f+D*beta block equivalent to the mixed R list matrix.
      integer, intent(out) :: info !! Zero on success; nonzero for size, syntax, coefficient, name, or start-vector errors.
      real(dp), intent(in), optional :: parameter_starts(:) !! Optional starts in first-occurrence parameter-name order.
      character(len=marss_parameter_name_len), allocatable :: names(:)
      real(dp), allocatable :: fixed(:)
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: starts(:)
      integer :: i
      integer :: parse_info

      info = 0
      if (size(numeric_values) /= size(expressions) .or. size(numeric_values) /= size(is_expression)) then
         info = 1
         return
      end if
      allocate(fixed(size(numeric_values)))
      fixed = numeric_values
      allocate(names(0))
      do i = 1, size(numeric_values)
         if (.not. is_expression(i)) cycle
         call scan_constraint_expression(expressions(i), fixed(i), names, parse_info)
         if (parse_info /= 0) then
            info = 10 + parse_info
            return
         end if
      end do
      allocate(design(size(numeric_values), size(names)))
      design = 0.0_dp
      do i = 1, size(numeric_values)
         if (.not. is_expression(i)) cycle
         call fill_constraint_expression(expressions(i), fixed(i), names, design(i, :), parse_info)
         if (parse_info /= 0) then
            info = 20 + parse_info
            return
         end if
      end do
      allocate(starts(size(names)))
      starts = 0.0_dp
      if (present(parameter_starts)) then
         if (size(parameter_starts) /= size(names)) then
            info = 2
            return
         end if
         starts = parameter_starts
      end if
      call marss_constraint_from_affine(fixed, design, starts, block, parse_info, names)
      if (parse_info /= 0) info = 30 + parse_info
   end subroutine marss_constraint_from_entries

   subroutine marss_constraint_from_labels(values, labels, block, info, parameter_names)
      real(dp), intent(in) :: values(:) !! Fixed values or numerical starting values for free labeled elements.
      integer, intent(in) :: labels(:) !! Zero marks fixed elements; positive contiguous labels identify shared free parameters.
      type(marss_constraint_block), intent(out) :: block !! Affine block equivalent to the numeric/string R list-matrix pattern.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid sizes, labels, or supplied parameter names.
      character(len=*), intent(in), optional :: parameter_names(:) !! Optional names for labels 1 through max(labels).
      integer, allocatable :: counts(:)
      integer :: i
      integer :: label
      integer :: p

      info = 0
      if (size(values) /= size(labels)) then
         info = 1
         return
      end if
      if (any(labels < 0)) then
         info = 2
         return
      end if
      p = 0
      if (size(labels) > 0) p = maxval(labels)
      if (present(parameter_names)) then
         if (size(parameter_names) /= p) then
            info = 3
            return
         end if
         do i = 1, p
            if (len_trim(parameter_names(i)) == 0) then
               info = 5
               return
            end if
            if (i > 1) then
               if (any(parameter_names(1:i - 1) == parameter_names(i))) then
                  info = 6
                  return
               end if
            end if
         end do
      end if
      allocate(block%fixed(size(values)), block%design(size(values), p), block%start(p))
      allocate(block%free_names(p))
      block%fixed = 0.0_dp
      block%design = 0.0_dp
      block%start = 0.0_dp
      allocate(counts(p))
      counts = 0
      do i = 1, size(values)
         label = labels(i)
         if (label == 0) then
            block%fixed(i) = values(i)
         else
            block%design(i, label) = 1.0_dp
            block%start(label) = block%start(label) + values(i)
            counts(label) = counts(label) + 1
         end if
      end do
      if (p > 0) then
         if (any(counts == 0)) then
            info = 4
            return
         end if
         do label = 1, p
            block%start(label) = block%start(label) / real(counts(label), dp)
            if (present(parameter_names)) then
               block%free_names(label) = trim(parameter_names(label))
            else
               write (block%free_names(label), '(a,i0)') 'beta', label
            end if
         end do
      end if
   end subroutine marss_constraint_from_labels

   subroutine marss_free_parameter_names(constraints, names)
      type(marss_constraints), intent(in) :: constraints !! Constraint blocks whose free-coordinate names are requested.
      character(len=marss_parameter_name_len), allocatable, intent(out) :: names(:) !! Names in standard beta-vector order.
      integer :: k

      allocate(names(marss_constraints_parameter_count(constraints)))
      k = 0
      call append_block_names(constraints%b, 'B', names, k)
      call append_block_names(constraints%u, 'U', names, k)
      call append_block_names(constraints%q, 'Q', names, k)
      call append_block_names(constraints%z, 'Z', names, k)
      call append_block_names(constraints%a, 'A', names, k)
      call append_block_names(constraints%r, 'R', names, k)
      call append_block_names(constraints%x0, 'x0', names, k)
      call append_block_names(constraints%v0, 'V0', names, k)
      call append_block_names(constraints%c, 'C', names, k)
      call append_block_names(constraints%d, 'D', names, k)
      call append_block_names(constraints%g, 'G', names, k)
      call append_block_names(constraints%h, 'H', names, k)
      call append_block_names(constraints%l, 'L', names, k)
   end subroutine marss_free_parameter_names

   subroutine marss_vectorized_parameter_names(constraints, names)
      type(marss_constraints), intent(in) :: constraints !! Constraint blocks whose MARSSvectorizeparam-style names are requested.
      character(len=marss_parameter_name_len), allocatable, intent(out) :: names(:) !! Block-prefixed names in beta-vector order.
      integer :: k

      allocate(names(marss_constraints_parameter_count(constraints)))
      k = 0
      call append_vectorized_block_names(constraints%b, 'B', names, k)
      call append_vectorized_block_names(constraints%u, 'U', names, k)
      call append_vectorized_block_names(constraints%q, 'Q', names, k)
      call append_vectorized_block_names(constraints%z, 'Z', names, k)
      call append_vectorized_block_names(constraints%a, 'A', names, k)
      call append_vectorized_block_names(constraints%r, 'R', names, k)
      call append_vectorized_block_names(constraints%x0, 'x0', names, k)
      call append_vectorized_block_names(constraints%v0, 'V0', names, k)
      call append_vectorized_block_names(constraints%c, 'C', names, k)
      call append_vectorized_block_names(constraints%d, 'D', names, k)
      call append_vectorized_block_names(constraints%g, 'G', names, k)
      call append_vectorized_block_names(constraints%h, 'H', names, k)
      call append_vectorized_block_names(constraints%l, 'L', names, k)
   end subroutine marss_vectorized_parameter_names

   subroutine marss_reorder_free_parameters(constraints, input_names, input_values, theta, info)
      type(marss_constraints), intent(in) :: constraints !! Constraints defining the canonical free-coordinate name order.
      character(len=*), intent(in) :: input_names(:) !! Raw or block-prefixed names associated with input_values.
      real(dp), intent(in) :: input_values(:) !! Free-parameter values in arbitrary named order.
      real(dp), allocatable, intent(out) :: theta(:) !! Values reordered to B,U,Q,Z,A,R,x0,V0,C,D,G,H,L beta order.
      integer, intent(out) :: info !! Zero on success; nonzero for size mismatch, duplicate, missing, or ambiguous names.
      character(len=marss_parameter_name_len), allocatable :: raw_names(:)
      character(len=marss_parameter_name_len), allocatable :: vector_names(:)
      logical, allocatable :: used(:)
      integer :: candidate
      integer :: i
      integer :: j
      integer :: match_count
      integer :: raw_count

      info = 0
      call marss_free_parameter_names(constraints, raw_names)
      call marss_vectorized_parameter_names(constraints, vector_names)
      if (size(input_names) /= size(vector_names) .or. size(input_values) /= size(vector_names)) then
         info = 1
         return
      end if
      allocate(theta(size(vector_names)), used(size(vector_names)))
      used = .false.
      do i = 1, size(vector_names)
         match_count = 0
         j = 0
         do candidate = 1, size(input_names)
            if (trim(input_names(candidate)) == trim(vector_names(i))) then
               match_count = match_count + 1
               j = candidate
            end if
         end do
         if (match_count == 0) then
            raw_count = count([(trim(raw_names(candidate)) == trim(raw_names(i)), candidate = 1, size(raw_names))])
            if (raw_count == 1) then
               do candidate = 1, size(input_names)
                  if (trim(input_names(candidate)) == trim(raw_names(i))) then
                     match_count = match_count + 1
                     j = candidate
                  end if
               end do
            end if
         end if
         if (match_count == 0) then
            info = 2
            return
         end if
         if (match_count > 1) then
            info = 3
            return
         end if
         if (used(j)) then
            info = 4
            return
         end if
         theta(i) = input_values(j)
         used(j) = .true.
      end do
      if (.not. all(used)) info = 5
   end subroutine marss_reorder_free_parameters

   subroutine marss_constraints_set_start_named(constraints, parameter_names, values, updated, info)
      type(marss_constraints), intent(in) :: constraints !! Constraint design and existing free-coordinate names.
      character(len=*), intent(in) :: parameter_names(:) !! Names associated with values in arbitrary input order.
      real(dp), intent(in) :: values(:) !! Starting beta values associated with parameter_names.
      type(marss_constraints), intent(out) :: updated !! Constraint copy with start coordinates replaced by named values.
      integer, intent(out) :: info !! Zero on success; otherwise a name-reordering or coordinate-count error.
      real(dp), allocatable :: theta(:)
      integer :: reorder_info

      call marss_reorder_free_parameters(constraints, parameter_names, values, theta, reorder_info)
      if (reorder_info /= 0) then
         info = reorder_info
         return
      end if
      call marss_constraints_set_start(constraints, theta, updated, info)
   end subroutine marss_constraints_set_start_named


   subroutine marss_constraints_update_start_named(constraints, parameter_names, values, updated, info)
      type(marss_constraints), intent(in) :: constraints !! Constraint design whose existing starts supply unspecified values.
      character(len=*), intent(in) :: parameter_names(:) !! Raw or block-prefixed names for the coordinates to override.
      real(dp), intent(in) :: values(:) !! Replacement beta starts corresponding one-to-one with parameter_names.
      type(marss_constraints), intent(out) :: updated !! Constraint copy with only the named start coordinates replaced.
      integer, intent(out) :: info !! Zero on success; nonzero for size mismatch, missing, duplicate, or ambiguous names.
      character(len=marss_parameter_name_len), allocatable :: raw_names(:)
      character(len=marss_parameter_name_len), allocatable :: vector_names(:)
      logical, allocatable :: changed(:)
      real(dp), allocatable :: theta(:)
      integer :: coord
      integer :: i
      integer :: j
      integer :: match_count
      integer :: raw_count

      info = 0
      if (size(parameter_names) /= size(values)) then
         info = 1
         return
      end if
      call marss_free_parameter_names(constraints, raw_names)
      call marss_vectorized_parameter_names(constraints, vector_names)
      call marss_constraint_start_vector(constraints, theta)
      allocate(changed(size(theta)))
      changed = .false.
      do i = 1, size(parameter_names)
         match_count = 0
         coord = 0
         do raw_count = 1, size(vector_names)
            if (trim(parameter_names(i)) == trim(vector_names(raw_count))) then
               match_count = match_count + 1
               coord = raw_count
            end if
         end do
         if (match_count == 0) then
            raw_count = count([(trim(raw_names(j)) == trim(parameter_names(i)), j = 1, size(raw_names))])
            if (raw_count == 1) then
               do j = 1, size(raw_names)
                  if (trim(parameter_names(i)) == trim(raw_names(j))) then
                     match_count = 1
                     coord = j
                     exit
                  end if
               end do
            else
               coord = 0
            end if
         end if
         if (match_count == 0 .or. coord == 0) then
            info = 2
            return
         end if
         if (match_count > 1) then
            info = 3
            return
         end if
         if (changed(coord)) then
            info = 4
            return
         end if
         theta(coord) = values(i)
         changed(coord) = .true.
      end do
      call marss_constraints_set_start(constraints, theta, updated, info)
   end subroutine marss_constraints_update_start_named

   pure subroutine marss_constraint_start_vector(constraints, theta)
      type(marss_constraints), intent(in) :: constraints !! Affine constraint blocks supplying starting beta vectors.
      real(dp), allocatable, intent(out) :: theta(:) !! Beta coordinates in B,U,Q,Z,A,R,x0,V0,C,D,G,H,L order.
      integer :: k

      allocate(theta(marss_constraints_parameter_count(constraints)))
      k = 0
      call append_block_start(constraints%b, theta, k)
      call append_block_start(constraints%u, theta, k)
      call append_block_start(constraints%q, theta, k)
      call append_block_start(constraints%z, theta, k)
      call append_block_start(constraints%a, theta, k)
      call append_block_start(constraints%r, theta, k)
      call append_block_start(constraints%x0, theta, k)
      call append_block_start(constraints%v0, theta, k)
      call append_block_start(constraints%c, theta, k)
      call append_block_start(constraints%d, theta, k)
      call append_block_start(constraints%g, theta, k)
      call append_block_start(constraints%h, theta, k)
      call append_block_start(constraints%l, theta, k)
   end subroutine marss_constraint_start_vector

   pure subroutine marss_constraints_set_start(constraints, theta, updated, info)
      type(marss_constraints), intent(in) :: constraints !! Constraint design and fixed blocks copied without structural changes.
      real(dp), intent(in) :: theta(:) !! Concatenated beta coordinates in B,U,Q,Z,A,R,x0,V0,C,D,G,H,L order.
      type(marss_constraints), intent(out) :: updated !! Copy with each specified block start vector replaced from theta.
      integer, intent(out) :: info !! Zero on success, or one when theta has the wrong number of free coordinates.
      integer :: k

      updated = constraints
      if (size(theta) /= marss_constraints_parameter_count(constraints)) then
         info = 1
         return
      end if
      k = 0
      call set_block_start(updated%b, theta, k)
      call set_block_start(updated%u, theta, k)
      call set_block_start(updated%q, theta, k)
      call set_block_start(updated%z, theta, k)
      call set_block_start(updated%a, theta, k)
      call set_block_start(updated%r, theta, k)
      call set_block_start(updated%x0, theta, k)
      call set_block_start(updated%v0, theta, k)
      call set_block_start(updated%c, theta, k)
      call set_block_start(updated%d, theta, k)
      call set_block_start(updated%g, theta, k)
      call set_block_start(updated%h, theta, k)
      call set_block_start(updated%l, theta, k)
      info = 0
   end subroutine marss_constraints_set_start

   pure subroutine marss_vectorize_free(model, constraints, theta, info)
      type(marss_model), intent(in) :: model !! Numerical model whose estimated affine beta coordinates are recovered.
      type(marss_constraints), intent(in) :: constraints !! Affine f+D*beta blocks defining fixed, free, and equality constraints.
      real(dp), allocatable, intent(out) :: theta(:) !! Concatenated beta vector in B,U,Q,Z,A,R,x0,V0,C,D,G,H,L order.
      integer, intent(out) :: info !! Zero on success; positive when values violate or do not identify the constraints.
      real(dp), allocatable :: values(:)
      integer :: k

      info = 0
      allocate(theta(marss_constraints_parameter_count(constraints)))
      k = 0

      if (allocated(model%b_t)) then
         values = reshape(model%b_t, [size(model%b_t)])
      else
         values = reshape(model%b, [size(model%b)])
      end if
      call project_block(constraints%b, values, theta, k, info)
      if (info /= 0) then
         info = 10 + info
         return
      end if

      if (allocated(model%u_t)) then
         values = reshape(model%u_t, [size(model%u_t)])
      else
         values = model%u
      end if
      call project_block(constraints%u, values, theta, k, info)
      if (info /= 0) then
         info = 20 + info
         return
      end if

      if (allocated(model%q_noise)) then
         if (allocated(model%q_noise_t)) then
            values = reshape(model%q_noise_t, [size(model%q_noise_t)])
         else
            values = reshape(model%q_noise, [size(model%q_noise)])
         end if
      else if (allocated(model%q_t)) then
         values = reshape(model%q_t, [size(model%q_t)])
      else
         values = reshape(model%q, [size(model%q)])
      end if
      call project_block(constraints%q, values, theta, k, info)
      if (info /= 0) then
         info = 30 + info
         return
      end if

      if (allocated(model%z_t)) then
         values = reshape(model%z_t, [size(model%z_t)])
      else
         values = reshape(model%z, [size(model%z)])
      end if
      call project_block(constraints%z, values, theta, k, info)
      if (info /= 0) then
         info = 40 + info
         return
      end if

      if (allocated(model%a_t)) then
         values = reshape(model%a_t, [size(model%a_t)])
      else
         values = model%a
      end if
      call project_block(constraints%a, values, theta, k, info)
      if (info /= 0) then
         info = 50 + info
         return
      end if

      if (allocated(model%r_noise)) then
         if (allocated(model%r_noise_t)) then
            values = reshape(model%r_noise_t, [size(model%r_noise_t)])
         else
            values = reshape(model%r_noise, [size(model%r_noise)])
         end if
      else if (allocated(model%r_t)) then
         values = reshape(model%r_t, [size(model%r_t)])
      else
         values = reshape(model%r, [size(model%r)])
      end if
      call project_block(constraints%r, values, theta, k, info)
      if (info /= 0) then
         info = 60 + info
         return
      end if

      values = model%x0
      call project_block(constraints%x0, values, theta, k, info)
      if (info /= 0) then
         info = 70 + info
         return
      end if

      if (allocated(model%v0_noise)) then
         values = reshape(model%v0_noise, [size(model%v0_noise)])
      else
         values = reshape(model%v0, [size(model%v0)])
      end if
      call project_block(constraints%v0, values, theta, k, info)
      if (info /= 0) then
         info = 80 + info
         return
      end if

      if (allocated(model%c_coef)) then
         if (allocated(model%c_coef_t)) then
            values = reshape(model%c_coef_t, [size(model%c_coef_t)])
         else
            values = reshape(model%c_coef, [size(model%c_coef)])
         end if
      else
         values = [real(dp) ::]
      end if
      call project_block(constraints%c, values, theta, k, info)
      if (info /= 0) then
         info = 90 + info
         return
      end if

      if (allocated(model%d_coef)) then
         if (allocated(model%d_coef_t)) then
            values = reshape(model%d_coef_t, [size(model%d_coef_t)])
         else
            values = reshape(model%d_coef, [size(model%d_coef)])
         end if
      else
         values = [real(dp) ::]
      end if
      call project_block(constraints%d, values, theta, k, info)
      if (info /= 0) then
         info = 100 + info
         return
      end if

      if (allocated(model%g)) then
         if (allocated(model%g_t)) then
            values = reshape(model%g_t, [size(model%g_t)])
         else
            values = reshape(model%g, [size(model%g)])
         end if
      else
         values = [real(dp) ::]
      end if
      call project_block(constraints%g, values, theta, k, info)
      if (info /= 0) then
         info = 110 + info
         return
      end if

      if (allocated(model%h)) then
         if (allocated(model%h_t)) then
            values = reshape(model%h_t, [size(model%h_t)])
         else
            values = reshape(model%h, [size(model%h)])
         end if
      else
         values = [real(dp) ::]
      end if
      call project_block(constraints%h, values, theta, k, info)
      if (info /= 0) then
         info = 120 + info
         return
      end if

      if (allocated(model%l)) then
         values = reshape(model%l, [size(model%l)])
      else
         values = [real(dp) ::]
      end if
      call project_block(constraints%l, values, theta, k, info)
      if (info /= 0) then
         info = 130 + info
         return
      end if
      if (k /= size(theta)) info = 140
   end subroutine marss_vectorize_free

   pure subroutine marss_unvectorize_free(template, constraints, theta, model, info)
      type(marss_model), intent(in) :: template !! Template carrying data, dimensions, and fixed numerical structure.
      type(marss_constraints), intent(in) :: constraints !! Affine f+D*beta blocks used to reconstruct model parameters.
      real(dp), intent(in) :: theta(:) !! Concatenated beta coordinates in the standard MARSS block order.
      type(marss_model), intent(out) :: model !! Reconstructed numerical model with the supplied free coordinates applied.
      integer, intent(out) :: info !! Zero on success; nonzero when dimensions or resulting covariance structure are invalid.

      call marss_apply_constraints(template, constraints, theta, model, info)
   end subroutine marss_unvectorize_free

   pure subroutine marss_apply_constraints(template, constraints, theta, model, info)
      type(marss_model), intent(in) :: template !! Template model whose unspecified blocks and observations are preserved.
      type(marss_constraints), intent(in) :: constraints !! Affine f+D*beta specifications for selected model blocks.
      real(dp), intent(in) :: theta(:) !! Concatenated beta vector matching marss_constraint_start_vector ordering.
      type(marss_model), intent(out) :: model !! Reconstructed numerical model after applying every specified affine block.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid coordinates, shapes, symmetry, or covariance.
      real(dp), allocatable :: values(:)
      integer :: k

      model = template
      if (size(theta) /= marss_constraints_parameter_count(constraints)) then
         info = 1
         return
      end if
      k = 0
      call block_values(constraints%b, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_b_size(model)) then
            info = 2
            return
         end if
         if (allocated(model%b_t)) then
            model%b_t = reshape(values, shape(model%b_t))
         else
            model%b = reshape(values, shape(model%b))
         end if
         deallocate(values)
      end if
      call block_values(constraints%u, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_u_size(model)) then
            info = 3
            return
         end if
         if (allocated(model%u_t)) then
            model%u_t = reshape(values, shape(model%u_t))
         else
            model%u = values
         end if
         deallocate(values)
      end if
      call block_values(constraints%q, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_q_size(model)) then
            info = 4
            return
         end if
         if (allocated(model%q_noise)) then
            if (allocated(model%q_noise_t)) then
               model%q_noise_t = reshape(values, shape(model%q_noise_t))
            else
               model%q_noise = reshape(values, shape(model%q_noise))
            end if
         else if (allocated(model%q_t)) then
            model%q_t = reshape(values, shape(model%q_t))
         else
            model%q = reshape(values, shape(model%q))
         end if
         deallocate(values)
      end if
      call block_values(constraints%z, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_z_size(model)) then
            info = 5
            return
         end if
         if (allocated(model%z_t)) then
            model%z_t = reshape(values, shape(model%z_t))
         else
            model%z = reshape(values, shape(model%z))
         end if
         deallocate(values)
      end if
      call block_values(constraints%a, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_a_size(model)) then
            info = 6
            return
         end if
         if (allocated(model%a_t)) then
            model%a_t = reshape(values, shape(model%a_t))
         else
            model%a = values
         end if
         deallocate(values)
      end if
      call block_values(constraints%r, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_r_size(model)) then
            info = 7
            return
         end if
         if (allocated(model%r_noise)) then
            if (allocated(model%r_noise_t)) then
               model%r_noise_t = reshape(values, shape(model%r_noise_t))
            else
               model%r_noise = reshape(values, shape(model%r_noise))
            end if
         else if (allocated(model%r_t)) then
            model%r_t = reshape(values, shape(model%r_t))
         else
            model%r = reshape(values, shape(model%r))
         end if
         deallocate(values)
      end if
      call block_values(constraints%x0, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= size(model%x0)) then
            info = 8
            return
         end if
         model%x0 = values
         deallocate(values)
      end if
      call block_values(constraints%v0, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_v0_size(model)) then
            info = 9
            return
         end if
         if (allocated(model%v0_noise)) then
            model%v0_noise = reshape(values, shape(model%v0_noise))
         else
            model%v0 = reshape(values, shape(model%v0))
         end if
         deallocate(values)
      end if
      call block_values(constraints%c, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_c_size(model)) then
            info = 10
            return
         end if
         if (allocated(model%c_coef_t)) then
            model%c_coef_t = reshape(values, shape(model%c_coef_t))
         else
            model%c_coef = reshape(values, shape(model%c_coef))
         end if
         deallocate(values)
      end if
      call block_values(constraints%d, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_d_size(model)) then
            info = 11
            return
         end if
         if (allocated(model%d_coef_t)) then
            model%d_coef_t = reshape(values, shape(model%d_coef_t))
         else
            model%d_coef = reshape(values, shape(model%d_coef))
         end if
         deallocate(values)
      end if
      call block_values(constraints%g, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_g_size(model)) then
            info = 12
            return
         end if
         if (allocated(model%g_t)) then
            model%g_t = reshape(values, shape(model%g_t))
         else
            model%g = reshape(values, shape(model%g))
         end if
         deallocate(values)
      end if
      call block_values(constraints%h, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_h_size(model)) then
            info = 13
            return
         end if
         if (allocated(model%h_t)) then
            model%h_t = reshape(values, shape(model%h_t))
         else
            model%h = reshape(values, shape(model%h))
         end if
         deallocate(values)
      end if
      call block_values(constraints%l, theta, k, values)
      if (allocated(values)) then
         if (size(values) /= active_l_size(model)) then
            info = 14
            return
         end if
         model%l = reshape(values, shape(model%l))
         deallocate(values)
      end if
      if (.not. marss_model_valid(model)) then
         info = 20
         return
      end if
      info = 0
   end subroutine marss_apply_constraints

   subroutine marss_optim_linear(start_model, constraints, max_iter, tol, fit)
      type(marss_model), intent(in) :: start_model !! Starting model; unspecified constraint blocks remain fixed at these values.
      type(marss_constraints), intent(in) :: constraints !! Vectorized affine f+D*beta constraints for estimated blocks.
      integer, intent(in) :: max_iter !! Maximum finite-difference BFGS iterations, which must be positive.
      real(dp), intent(in) :: tol !! Relative objective and gradient convergence tolerance, which must be positive.
      type(marss_fit_result), intent(out) :: fit !! Constrained optimum, Kalman result, convergence flag, and log likelihood.
      type(marss_model) :: final_model
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: theta_new(:)
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: gradient_new(:)
      real(dp), allocatable :: h_inv(:, :)
      real(dp), allocatable :: direction(:)
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: left(:, :)
      real(dp), allocatable :: right(:, :)
      real(dp) :: alpha
      real(dp) :: f
      real(dp) :: f_new
      real(dp) :: gtd
      real(dp) :: rho
      real(dp) :: ys
      integer :: info
      integer :: iter
      integer :: p

      fit%converged = .false.
      fit%iterations = 0
      fit%info = 0
      if (max_iter < 1 .or. tol <= 0.0_dp) then
         fit%info = 1
         return
      end if
      if (.not. marss_constraints_valid(start_model, constraints)) then
         fit%info = 2
         return
      end if
      call marss_constraint_start_vector(constraints, theta)
      p = size(theta)
      if (p == 0) then
         fit%model = start_model
         call marss_kf(fit%model, fit%kf)
         if (.not. fit%kf%ok) then
            fit%info = 20 + fit%kf%info
            return
         end if
         fit%loglik = fit%kf%loglik
         allocate(fit%free_parameters(0))
         fit%converged = .true.
         return
      end if
      allocate(theta_new(p), gradient(p), gradient_new(p), h_inv(p, p), direction(p), s(p), y(p))
      allocate(left(p, p), right(p, p))
      call constrained_objective(start_model, constraints, theta, f, info)
      if (info /= 0) then
         fit%info = 30 + info
         return
      end if
      call constrained_gradient(start_model, constraints, theta, f, gradient, info)
      if (info /= 0) then
         fit%info = 40 + info
         return
      end if
      h_inv = identity_matrix(p)

      do iter = 1, max_iter
         if (maxval(abs(gradient)) <= tol) then
            fit%converged = .true.
            exit
         end if
         direction = -matmul(h_inv, gradient)
         gtd = dot_product(gradient, direction)
         if (gtd >= -epsilon(1.0_dp)) then
            h_inv = identity_matrix(p)
            direction = -gradient
         end if
         call constrained_armijo(start_model, constraints, theta, f, gradient, direction, theta_new, f_new, alpha, info)
         if (info /= 0) then
            fit%info = 50 + info
            return
         end if
         call constrained_gradient(start_model, constraints, theta_new, f_new, gradient_new, info)
         if (info /= 0) then
            fit%info = 60 + info
            return
         end if
         s = theta_new - theta
         y = gradient_new - gradient
         ys = dot_product(y, s)
         if (ys > sqrt(epsilon(1.0_dp)) * max(1.0_dp, sqrt(dot_product(y, y) * dot_product(s, s)))) then
            rho = 1.0_dp / ys
            left = identity_matrix(p) - rho * outer_product(s, y)
            right = identity_matrix(p) - rho * outer_product(y, s)
            h_inv = matmul(matmul(left, h_inv), right) + rho * outer_product(s, s)
         else
            h_inv = identity_matrix(p)
         end if
         if (abs(f_new - f) <= tol * (1.0_dp + abs(f))) fit%converged = .true.
         theta = theta_new
         gradient = gradient_new
         f = f_new
         if (fit%converged) exit
      end do
      fit%iterations = min(iter, max_iter)
      call marss_apply_constraints(start_model, constraints, theta, final_model, info)
      if (info /= 0) then
         fit%info = 70 + info
         return
      end if
      fit%model = final_model
      call marss_kf(fit%model, fit%kf)
      if (.not. fit%kf%ok) then
         fit%info = 80 + fit%kf%info
         return
      end if
      fit%loglik = fit%kf%loglik
      fit%free_parameters = theta
      fit%info = 0
   end subroutine marss_optim_linear

   subroutine marss_hessian_linear(template, constraints, theta, hessian, rel_step, info)
      type(marss_model), intent(in) :: template !! Template model defining data, fixed blocks, and active constrained shapes.
      type(marss_constraints), intent(in) :: constraints !! Affine f+D*beta constraints defining the free-coordinate model.
      real(dp), intent(in) :: theta(:) !! Free beta coordinates at which the log-likelihood Hessian is evaluated.
      real(dp), allocatable, intent(out) :: hessian(:, :) !! Central finite-difference Hessian with respect to beta.
      real(dp), intent(in), optional :: rel_step !! Relative finite-difference step, defaulting to 1e-4.
      integer, intent(out) :: info !! Zero on success, nonzero if the base or a perturbed constrained model is invalid.
      real(dp), allocatable :: work(:)
      real(dp) :: f0
      real(dp) :: fm
      real(dp) :: fmm
      real(dp) :: fmp
      real(dp) :: fp
      real(dp) :: fpm
      real(dp) :: fpp
      real(dp) :: hi
      real(dp) :: hj
      real(dp) :: step
      integer :: i
      integer :: j
      integer :: p

      p = size(theta)
      if (p /= marss_constraints_parameter_count(constraints)) then
         allocate(hessian(0, 0))
         info = 1
         return
      end if
      allocate(hessian(p, p), work(p))
      if (p == 0) then
         info = 0
         return
      end if
      step = 1.0e-4_dp
      if (present(rel_step)) step = rel_step
      if (step <= 0.0_dp) then
         info = 2
         return
      end if
      call constrained_loglik(template, constraints, theta, f0, info)
      if (info /= 0) return
      hessian = 0.0_dp
      do i = 1, p
         hi = step * max(1.0_dp, abs(theta(i)))
         work = theta
         work(i) = theta(i) + hi
         call constrained_loglik(template, constraints, work, fp, info)
         if (info /= 0) return
         work(i) = theta(i) - hi
         call constrained_loglik(template, constraints, work, fm, info)
         if (info /= 0) return
         hessian(i, i) = (fp - 2.0_dp*f0 + fm) / (hi*hi)
         do j = i + 1, p
            hj = step * max(1.0_dp, abs(theta(j)))
            work = theta
            work(i) = theta(i) + hi
            work(j) = theta(j) + hj
            call constrained_loglik(template, constraints, work, fpp, info)
            if (info /= 0) return
            work(j) = theta(j) - hj
            call constrained_loglik(template, constraints, work, fpm, info)
            if (info /= 0) return
            work(i) = theta(i) - hi
            work(j) = theta(j) + hj
            call constrained_loglik(template, constraints, work, fmp, info)
            if (info /= 0) return
            work(j) = theta(j) - hj
            call constrained_loglik(template, constraints, work, fmm, info)
            if (info /= 0) return
            hessian(i, j) = (fpp - fpm - fmp + fmm) / (4.0_dp*hi*hj)
            hessian(j, i) = hessian(i, j)
         end do
      end do
      info = 0
   end subroutine marss_hessian_linear

   subroutine marss_fisher_i_linear(template, constraints, theta, fisher, rel_step, info)
      type(marss_model), intent(in) :: template !! Template model defining the constrained likelihood.
      type(marss_constraints), intent(in) :: constraints !! Affine constraints defining free beta coordinates.
      real(dp), intent(in) :: theta(:) !! Free beta coordinates at which observed information is evaluated.
      real(dp), allocatable, intent(out) :: fisher(:, :) !! Negative beta-coordinate log-likelihood Hessian.
      real(dp), intent(in), optional :: rel_step !! Relative finite-difference step, defaulting to 1e-4.
      integer, intent(out) :: info !! Zero on success, nonzero when Hessian evaluation fails.
      real(dp), allocatable :: hessian(:, :)

      if (present(rel_step)) then
         call marss_hessian_linear(template, constraints, theta, hessian, rel_step, info)
      else
         call marss_hessian_linear(template, constraints, theta, hessian, info=info)
      end if
      if (info == 0) then
         allocate(fisher(size(hessian, 1), size(hessian, 2)))
         fisher = -hessian
      else
         allocate(fisher(0, 0))
      end if
   end subroutine marss_fisher_i_linear

   subroutine marss_param_cis_linear(template, constraints, theta, alpha, se, lower, upper, info, rel_step)
      type(marss_model), intent(in) :: template !! Template model defining the constrained likelihood.
      type(marss_constraints), intent(in) :: constraints !! Affine constraints defining free beta coordinates.
      real(dp), intent(in) :: theta(:) !! Estimated free beta coordinates.
      real(dp), intent(in) :: alpha !! Two-sided significance level, for example 0.05 for 95 percent intervals.
      real(dp), allocatable, intent(out) :: se(:) !! Hessian-based standard errors for beta.
      real(dp), allocatable, intent(out) :: lower(:) !! Lower normal-approximation beta confidence limits.
      real(dp), allocatable, intent(out) :: upper(:) !! Upper normal-approximation beta confidence limits.
      integer, intent(out) :: info !! Zero on success, nonzero if observed information is singular or invalid.
      real(dp), intent(in), optional :: rel_step !! Relative finite-difference step, defaulting to 1e-4.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fisher(:, :)
      real(dp) :: zcrit
      integer :: i

      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         info = 1
         return
      end if
      if (present(rel_step)) then
         call marss_fisher_i_linear(template, constraints, theta, fisher, rel_step, info)
      else
         call marss_fisher_i_linear(template, constraints, theta, fisher, info=info)
      end if
      if (info /= 0) return
      if (size(theta) == 0) then
         allocate(se(0), lower(0), upper(0))
         return
      end if
      call inverse_matrix(fisher, covariance, info)
      if (info /= 0) return
      allocate(se(size(theta)), lower(size(theta)), upper(size(theta)))
      do i = 1, size(theta)
         se(i) = sqrt(max(covariance(i, i), 0.0_dp))
      end do
      zcrit = normal_quantile(1.0_dp - 0.5_dp*alpha)
      lower = theta - zcrit*se
      upper = theta + zcrit*se
      info = 0
   end subroutine marss_param_cis_linear

   subroutine constrained_loglik(template, constraints, theta, value, info)
      type(marss_model), intent(in) :: template !! Template model used to reconstruct the constrained numerical model.
      type(marss_constraints), intent(in) :: constraints !! Affine constraints defining the beta coordinates.
      real(dp), intent(in) :: theta(:) !! Free beta coordinates at which log likelihood is evaluated.
      real(dp), intent(out) :: value !! Kalman innovations log likelihood.
      integer, intent(out) :: info !! Zero on success, nonzero when reconstruction or filtering fails.
      type(marss_model) :: model
      type(marss_kf_result) :: kf

      call marss_apply_constraints(template, constraints, theta, model, info)
      if (info /= 0) return
      call marss_kf(model, kf, smoother=.false.)
      if (.not. kf%ok) then
         info = 100 + kf%info
         return
      end if
      value = kf%loglik
      info = 0
   end subroutine constrained_loglik

   subroutine constrained_objective(template, constraints, theta, value, info)
      type(marss_model), intent(in) :: template !! Template model for reconstructing constrained optimizer coordinates.
      type(marss_constraints), intent(in) :: constraints !! Affine block constraints defining theta.
      real(dp), intent(in) :: theta(:) !! Constrained optimizer beta coordinates.
      real(dp), intent(out) :: value !! Negative Kalman innovations log likelihood, or huge on invalid covariance/model values.
      integer, intent(out) :: info !! Zero on success, nonzero if reconstruction or filtering fails.
      type(marss_model) :: model
      type(marss_kf_result) :: kf

      call marss_apply_constraints(template, constraints, theta, model, info)
      if (info /= 0) then
         value = huge(1.0_dp)
         return
      end if
      call marss_kf(model, kf, smoother=.false.)
      if (.not. kf%ok) then
         info = 100 + kf%info
         value = huge(1.0_dp)
         return
      end if
      value = -kf%loglik
      info = 0
   end subroutine constrained_objective

   subroutine constrained_gradient(template, constraints, theta, f0, gradient, info)
      type(marss_model), intent(in) :: template !! Template model for constrained finite-difference evaluations.
      type(marss_constraints), intent(in) :: constraints !! Affine block constraints defining theta.
      real(dp), intent(in) :: theta(:) !! Current constrained optimizer beta coordinates.
      real(dp), intent(in) :: f0 !! Objective at theta, used for one-sided fallback differences.
      real(dp), intent(out) :: gradient(:) !! Finite-difference gradient of the constrained negative log likelihood.
      integer, intent(out) :: info !! Zero on success, nonzero when both perturbation directions are invalid.
      real(dp), allocatable :: work(:)
      real(dp) :: fm
      real(dp) :: fp
      real(dp) :: h
      integer :: im
      integer :: ip
      integer :: j

      if (size(gradient) /= size(theta)) then
         info = 1
         return
      end if
      allocate(work(size(theta)))
      do j = 1, size(theta)
         h = 1.0e-5_dp * (1.0_dp + abs(theta(j)))
         work = theta
         work(j) = theta(j) + h
         call constrained_objective(template, constraints, work, fp, ip)
         work(j) = theta(j) - h
         call constrained_objective(template, constraints, work, fm, im)
         if (ip == 0 .and. im == 0) then
            gradient(j) = (fp - fm) / (2.0_dp * h)
         else if (ip == 0) then
            gradient(j) = (fp - f0) / h
         else if (im == 0) then
            gradient(j) = (f0 - fm) / h
         else
            info = 10 + j
            return
         end if
      end do
      info = 0
   end subroutine constrained_gradient

   subroutine constrained_armijo(template, constraints, theta, f, gradient, direction, theta_new, f_new, alpha, info)
      type(marss_model), intent(in) :: template !! Template model for constrained line-search objective evaluations.
      type(marss_constraints), intent(in) :: constraints !! Affine block constraints defining theta.
      real(dp), intent(in) :: theta(:) !! Current constrained optimizer coordinates.
      real(dp), intent(in) :: f !! Current negative log likelihood.
      real(dp), intent(in) :: gradient(:) !! Current objective gradient.
      real(dp), intent(in) :: direction(:) !! Descent direction proposed by constrained BFGS.
      real(dp), intent(out) :: theta_new(:) !! Accepted constrained optimizer coordinates.
      real(dp), intent(out) :: f_new !! Objective value at theta_new.
      real(dp), intent(out) :: alpha !! Accepted line-search step length.
      integer, intent(out) :: info !! Zero on success, nonzero if no covariance-valid Armijo step is found.
      real(dp), parameter :: c1 = 1.0e-4_dp
      real(dp) :: slope
      integer :: eval_info
      integer :: trial

      slope = dot_product(gradient, direction)
      alpha = 1.0_dp
      do trial = 1, 40
         theta_new = theta + alpha * direction
         call constrained_objective(template, constraints, theta_new, f_new, eval_info)
         if (eval_info == 0) then
            if (f_new <= f + c1 * alpha * slope) then
               info = 0
               return
            end if
         end if
         alpha = 0.5_dp * alpha
      end do
      theta_new = theta
      f_new = f
      info = 1
   end subroutine constrained_armijo

   pure subroutine project_block(block, values, theta, k, info)
      type(marss_constraint_block), intent(in) :: block !! Affine block whose beta coordinates are projected from numerical values.
      real(dp), intent(in) :: values(:) !! Flattened active numerical block in Fortran column-major order.
      real(dp), intent(inout) :: theta(:) !! Destination concatenated beta vector.
      integer, intent(inout) :: k !! Number of beta coordinates already filled and updated on return.
      integer, intent(out) :: info !! Zero on success; nonzero for shape, rank, or constraint-residual failures.
      real(dp), allocatable :: gram(:, :)
      real(dp), allocatable :: gram_inverse(:, :)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: residual(:)
      real(dp) :: tolerance
      integer :: inv_info
      integer :: p

      info = 0
      if (.not. allocated(block%design)) return
      if (.not. allocated(block%fixed)) then
         info = 1
         return
      end if
      if (size(values) /= size(block%fixed)) then
         info = 2
         return
      end if
      p = size(block%design, 2)
      tolerance = 1.0e-8_dp * (1.0_dp + max(0.0_dp, maxval(abs(values))))
      if (p == 0) then
         if (size(values) > 0) then
            if (maxval(abs(values - block%fixed)) > tolerance) info = 3
         end if
         return
      end if
      gram = matmul(transpose(block%design), block%design)
      call inverse_matrix(gram, gram_inverse, inv_info)
      if (inv_info /= 0) then
         info = 4
         return
      end if
      beta = matmul(gram_inverse, matmul(transpose(block%design), values - block%fixed))
      residual = block%fixed + matmul(block%design, beta) - values
      if (maxval(abs(residual)) > tolerance) then
         info = 5
         return
      end if
      theta(k + 1:k + p) = beta
      k = k + p
   end subroutine project_block

   pure logical function block_valid(block, nelem) result(ok)
      type(marss_constraint_block), intent(in) :: block !! Affine block checked for f, D, and beta dimension consistency.
      integer, intent(in) :: nelem !! Number of active numerical elements in the corresponding model block.

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
      if (allocated(block%free_names)) then
         if (size(block%free_names) /= size(block%start)) ok = .false.
      end if
   end function block_valid

   pure integer function block_parameter_count(block) result(nparam)
      type(marss_constraint_block), intent(in) :: block !! Affine block whose beta-coordinate count is requested.

      nparam = 0
      if (allocated(block%design)) nparam = size(block%design, 2)
   end function block_parameter_count

   pure subroutine set_block_start(block, theta, k)
      type(marss_constraint_block), intent(inout) :: block !! Constraint block whose stored starting beta coordinates are replaced.
      real(dp), intent(in) :: theta(:) !! Concatenated beta vector supplying this block's replacement start coordinates.
      integer, intent(inout) :: k !! Number of beta coordinates already consumed and updated on return.
      integer :: p

      if (.not. allocated(block%design)) return
      p = size(block%design, 2)
      block%start = theta(k + 1:k + p)
      k = k + p
   end subroutine set_block_start

   subroutine scan_constraint_expression(expression, fixed_value, names, info)
      character(len=*), intent(in) :: expression !! R-style character entry to scan for its fixed offset and parameter labels.
      real(dp), intent(out) :: fixed_value !! Parsed fixed offset, zero when the expression starts with a free term.
      character(len=marss_parameter_name_len), allocatable, intent(inout) :: names(:) !! Unique labels in first-occurrence order.
      integer, intent(out) :: info !! Zero on success; nonzero for malformed syntax, coefficient, or parameter name.
      character(len=:), allocatable :: text
      character(len=:), allocatable :: term
      integer :: first
      integer :: last
      integer :: next_plus
      logical :: first_term
      logical :: fixed_consumed

      info = 0
      fixed_value = 0.0_dp
      text = trim(adjustl(expression))
      if (len(text) == 0) then
         info = 1
         return
      end if
      if (index(text, '+') == 0 .and. index(text, '*') == 0) then
         call add_constraint_name(text, names, info)
         return
      end if
      first = 1
      first_term = .true.
      fixed_consumed = .false.
      do
         next_plus = index(text(first:), '+')
         if (next_plus == 0) then
            last = len(text)
         else
            last = first + next_plus - 2
         end if
         if (last < first) then
            info = 2
            return
         end if
         term = trim(adjustl(text(first:last)))
         if (len(term) == 0) then
            info = 2
            return
         end if
         call scan_constraint_term(term, first_term, fixed_consumed, fixed_value, names, info)
         if (info /= 0) return
         first_term = .false.
         if (next_plus == 0) exit
         first = last + 2
         if (first > len(text)) then
            info = 2
            return
         end if
      end do
   end subroutine scan_constraint_expression

   subroutine fill_constraint_expression(expression, fixed_value, names, row, info)
      character(len=*), intent(in) :: expression !! Validated R-style character entry to convert to one design-matrix row.
      real(dp), intent(in) :: fixed_value !! Fixed offset obtained during the scan pass; retained for interface clarity.
      character(len=marss_parameter_name_len), intent(in) :: names(:) !! Canonical parameter-label order for design columns.
      real(dp), intent(out) :: row(:) !! Design coefficients for this flattened model-matrix entry.
      integer, intent(out) :: info !! Zero on success; nonzero if syntax or a parameter label cannot be resolved.
      character(len=:), allocatable :: text
      character(len=:), allocatable :: term
      integer :: first
      integer :: last
      integer :: next_plus
      logical :: first_term
      logical :: fixed_consumed
      real(dp) :: ignored_fixed

      info = 0
      row = 0.0_dp
      ignored_fixed = fixed_value
      text = trim(adjustl(expression))
      if (index(text, '+') == 0 .and. index(text, '*') == 0) then
         call add_constraint_coefficient(text, 1.0_dp, names, row, info)
         return
      end if
      first = 1
      first_term = .true.
      fixed_consumed = .false.
      do
         next_plus = index(text(first:), '+')
         if (next_plus == 0) then
            last = len(text)
         else
            last = first + next_plus - 2
         end if
         term = trim(adjustl(text(first:last)))
         call fill_constraint_term(term, first_term, fixed_consumed, names, row, info)
         if (info /= 0) return
         first_term = .false.
         if (next_plus == 0) exit
         first = last + 2
      end do
   end subroutine fill_constraint_expression

   subroutine scan_constraint_term(term, first_term, fixed_consumed, fixed_value, names, info)
      character(len=*), intent(in) :: term !! Single plus-delimited term from an affine character expression.
      logical, intent(in) :: first_term !! True only for the first term, where a numeric fixed offset is permitted.
      logical, intent(inout) :: fixed_consumed !! Tracks whether the optional leading fixed offset was recognized.
      real(dp), intent(inout) :: fixed_value !! Parsed fixed offset for the expression.
      character(len=marss_parameter_name_len), allocatable, intent(inout) :: names(:) !! Unique free labels collected so far.
      integer, intent(out) :: info !! Zero on success; nonzero for malformed numeric coefficients or labels.
      character(len=:), allocatable :: name
      real(dp) :: coefficient
      real(dp) :: numeric_term
      integer :: ios
      integer :: star

      info = 0
      star = index(term, '*')
      if (first_term .and. star == 0) then
         read (term, *, iostat=ios) numeric_term
         if (ios == 0) then
            fixed_value = numeric_term
            fixed_consumed = .true.
            return
         end if
      end if
      if (star == 0) then
         if (.not. first_term) then
            read (term, *, iostat=ios) numeric_term
            if (ios == 0) then
               info = 8
               return
            end if
         end if
         coefficient = 1.0_dp
         name = trim(adjustl(term))
      else
         if (index(term(star + 1:), '*') /= 0 .or. star == 1 .or. star == len_trim(term)) then
            info = 3
            return
         end if
         read (term(:star - 1), *, iostat=ios) coefficient
         if (ios /= 0) then
            info = 4
            return
         end if
         name = trim(adjustl(term(star + 1:)))
      end if
      call add_constraint_name(name, names, info)
   end subroutine scan_constraint_term

   subroutine fill_constraint_term(term, first_term, fixed_consumed, names, row, info)
      character(len=*), intent(in) :: term !! Single plus-delimited term from a previously validated affine expression.
      logical, intent(in) :: first_term !! True for the first term, where a numeric fixed offset may be skipped.
      logical, intent(inout) :: fixed_consumed !! Tracks whether a leading numeric fixed term was skipped.
      character(len=marss_parameter_name_len), intent(in) :: names(:) !! Canonical parameter labels for the design row.
      real(dp), intent(inout) :: row(:) !! Design row receiving this term's coefficient, accumulated for repeated labels.
      integer, intent(out) :: info !! Zero on success; nonzero for coefficient syntax or an unresolved label.
      character(len=:), allocatable :: name
      real(dp) :: coefficient
      real(dp) :: numeric_term
      integer :: ios
      integer :: star

      info = 0
      star = index(term, '*')
      if (first_term .and. star == 0) then
         read (term, *, iostat=ios) numeric_term
         if (ios == 0) then
            fixed_consumed = .true.
            return
         end if
      end if
      if (star == 0) then
         if (.not. first_term) then
            read (term, *, iostat=ios) numeric_term
            if (ios == 0) then
               info = 8
               return
            end if
         end if
         coefficient = 1.0_dp
         name = trim(adjustl(term))
      else
         read (term(:star - 1), *, iostat=ios) coefficient
         if (ios /= 0) then
            info = 4
            return
         end if
         name = trim(adjustl(term(star + 1:)))
      end if
      call add_constraint_coefficient(name, coefficient, names, row, info)
   end subroutine fill_constraint_term

   subroutine add_constraint_name(name, names, info)
      character(len=*), intent(in) :: name !! Candidate free-parameter label from an R-style list-matrix entry.
      character(len=marss_parameter_name_len), allocatable, intent(inout) :: names(:) !! Unique labels in first-use order.
      integer, intent(out) :: info !! Zero on success; nonzero when the label is blank or exceeds the supported name length.
      character(len=marss_parameter_name_len), allocatable :: grown(:)
      integer :: i

      info = 0
      if (len_trim(name) == 0) then
         info = 5
         return
      end if
      if (len_trim(name) > marss_parameter_name_len) then
         info = 6
         return
      end if
      do i = 1, size(names)
         if (trim(names(i)) == trim(name)) return
      end do
      allocate(grown(size(names) + 1))
      if (size(names) > 0) grown(1:size(names)) = names
      grown(size(grown)) = trim(name)
      call move_alloc(grown, names)
   end subroutine add_constraint_name

   subroutine add_constraint_coefficient(name, coefficient, names, row, info)
      character(len=*), intent(in) :: name !! Parameter label whose coefficient is accumulated into the design row.
      real(dp), intent(in) :: coefficient !! Numeric multiplier on this occurrence of the free parameter.
      character(len=marss_parameter_name_len), intent(in) :: names(:) !! Canonical parameter labels for design columns.
      real(dp), intent(inout) :: row(:) !! Design row updated at the matching parameter column.
      integer, intent(out) :: info !! Zero on success; nonzero when the parameter label is not in the scanned name set.
      integer :: i

      info = 7
      do i = 1, size(names)
         if (trim(names(i)) == trim(name)) then
            row(i) = row(i) + coefficient
            info = 0
            return
         end if
      end do
   end subroutine add_constraint_coefficient

   subroutine append_vectorized_block_names(block, prefix, names, k)
      type(marss_constraint_block), intent(in) :: block !! One parameter block whose names receive the MARSS block prefix.
      character(len=*), intent(in) :: prefix !! MARSS parameter-block name such as B, Q, x0, or V0.
      character(len=marss_parameter_name_len), intent(inout) :: names(:) !! Destination vector-name array.
      integer, intent(inout) :: k !! Number of names already written before this block.
      character(len=marss_parameter_name_len) :: raw_name
      integer :: i
      integer :: p

      p = block_parameter_count(block)
      do i = 1, p
         k = k + 1
         raw_name = ''
         if (allocated(block%free_names)) then
            raw_name = trim(block%free_names(i))
         else
            write (raw_name, '(i0)') i
         end if
         names(k) = trim(prefix) // '.' // trim(raw_name)
      end do
   end subroutine append_vectorized_block_names

   subroutine append_block_names(block, prefix, names, k)
      type(marss_constraint_block), intent(in) :: block !! Affine block whose free-coordinate names are appended when specified.
      character(len=*), intent(in) :: prefix !! Fallback block prefix used when explicit names are unavailable.
      character(len=marss_parameter_name_len), intent(inout) :: names(:) !! Destination names in standard beta-vector order.
      integer, intent(inout) :: k !! Number of names already filled and updated on return.
      integer :: i
      integer :: p

      if (.not. allocated(block%design)) return
      p = size(block%design, 2)
      do i = 1, p
         if (allocated(block%free_names)) then
            if (len_trim(block%free_names(i)) > 0) then
               names(k + i) = trim(block%free_names(i))
            else
               write (names(k + i), '(a,i0)') trim(prefix), i
            end if
         else
            write (names(k + i), '(a,i0)') trim(prefix), i
         end if
      end do
      k = k + p
   end subroutine append_block_names

   pure subroutine append_block_start(block, theta, k)
      type(marss_constraint_block), intent(in) :: block !! Affine block whose starting beta vector is appended when specified.
      real(dp), intent(inout) :: theta(:) !! Destination concatenated beta vector.
      integer, intent(inout) :: k !! Number of beta coordinates already filled and updated on return.

      if (.not. allocated(block%design)) return
      theta(k + 1:k + size(block%start)) = block%start
      k = k + size(block%start)
   end subroutine append_block_start

   pure subroutine block_values(block, theta, k, values)
      type(marss_constraint_block), intent(in) :: block !! Affine f+D*beta block to evaluate when specified.
      real(dp), intent(in) :: theta(:) !! Concatenated beta vector supplying this block's coordinates.
      integer, intent(inout) :: k !! Number of beta coordinates already consumed and updated on return.
      real(dp), allocatable, intent(out) :: values(:) !! Evaluated flattened block, or unallocated when the block is unspecified.
      integer :: p

      if (.not. allocated(block%design)) return
      p = size(block%design, 2)
      allocate(values(size(block%fixed)))
      values = block%fixed + matmul(block%design, theta(k + 1:k + p))
      k = k + p
   end subroutine block_values

   pure integer function active_b_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active B representation size is requested.

      if (allocated(model%b_t)) then
         n = size(model%b_t)
      else
         n = size(model%b)
      end if
   end function active_b_size

   pure integer function active_u_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active U representation size is requested.

      if (allocated(model%u_t)) then
         n = size(model%u_t)
      else
         n = size(model%u)
      end if
   end function active_u_size

   pure integer function active_q_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active raw Q representation size is requested.

      if (allocated(model%q_noise)) then
         if (allocated(model%q_noise_t)) then
            n = size(model%q_noise_t)
         else
            n = size(model%q_noise)
         end if
      else if (allocated(model%q_t)) then
         n = size(model%q_t)
      else
         n = size(model%q)
      end if
   end function active_q_size

   pure integer function active_z_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active Z representation size is requested.

      if (allocated(model%z_t)) then
         n = size(model%z_t)
      else
         n = size(model%z)
      end if
   end function active_z_size

   pure integer function active_a_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active A representation size is requested.

      if (allocated(model%a_t)) then
         n = size(model%a_t)
      else
         n = size(model%a)
      end if
   end function active_a_size

   pure integer function active_r_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active raw R representation size is requested.

      if (allocated(model%r_noise)) then
         if (allocated(model%r_noise_t)) then
            n = size(model%r_noise_t)
         else
            n = size(model%r_noise)
         end if
      else if (allocated(model%r_t)) then
         n = size(model%r_t)
      else
         n = size(model%r)
      end if
   end function active_r_size

   pure integer function active_v0_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active raw V0 representation size is requested.

      if (allocated(model%v0_noise)) then
         n = size(model%v0_noise)
      else
         n = size(model%v0)
      end if
   end function active_v0_size

   pure integer function active_c_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active C coefficient representation size is requested.

      n = 0
      if (allocated(model%c_coef_t)) then
         n = size(model%c_coef_t)
      else if (allocated(model%c_coef)) then
         n = size(model%c_coef)
      end if
   end function active_c_size

   pure integer function active_d_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active D coefficient representation size is requested.

      n = 0
      if (allocated(model%d_coef_t)) then
         n = size(model%d_coef_t)
      else if (allocated(model%d_coef)) then
         n = size(model%d_coef)
      end if
   end function active_d_size

   pure integer function active_g_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active G loading representation size is requested.

      n = 0
      if (allocated(model%g_t)) then
         n = size(model%g_t)
      else if (allocated(model%g)) then
         n = size(model%g)
      end if
   end function active_g_size

   pure integer function active_h_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active H loading representation size is requested.

      n = 0
      if (allocated(model%h_t)) then
         n = size(model%h_t)
      else if (allocated(model%h)) then
         n = size(model%h)
      end if
   end function active_h_size

   pure integer function active_l_size(model) result(n)
      type(marss_model), intent(in) :: model !! Model whose active L loading representation size is requested.

      n = 0
      if (allocated(model%l)) n = size(model%l)
   end function active_l_size

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:) !! Left vector in the rank-one product.
      real(dp), intent(in) :: y(:) !! Right vector in the rank-one product.
      real(dp) :: a(size(x), size(y))

      a = spread(x, 2, size(y)) * spread(y, 1, size(x))
   end function outer_product

end module marss_constraints_mod
