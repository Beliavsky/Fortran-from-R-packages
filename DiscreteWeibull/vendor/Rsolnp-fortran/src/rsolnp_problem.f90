! SPDX-License-Identifier: GPL-2.0-only
module rsolnp_problem
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rsolnp_kinds, only : dp
   use rsolnp_types, only : solnp_problem, solnp_success, solnp_invalid_problem
   implicit none
   private

   real(dp), parameter :: default_bound = 1.0e20_dp

   public :: prepare_problem, solnp_standardize_problem

contains

   subroutine prepare_problem(problem, status, message)
      type(solnp_problem), intent(inout) :: problem
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      integer :: n

      status = solnp_invalid_problem
      message = ''
      if (.not. associated(problem%fn)) then
         message = 'objective callback is not associated'
         return
      end if
      if (.not. allocated(problem%start)) then
         message = 'starting parameters are not allocated'
         return
      end if
      n = size(problem%start)
      if (n < 1) then
         message = 'at least one parameter is required'
         return
      end if
      if (problem%n == 0) problem%n = n
      if (problem%n /= n) then
         message = 'problem%n differs from size(problem%start)'
         return
      end if

      if (.not. allocated(problem%lower)) then
         allocate(problem%lower(n))
         problem%lower = -default_bound
      end if
      if (.not. allocated(problem%upper)) then
         allocate(problem%upper(n))
         problem%upper = default_bound
      end if
      if (size(problem%lower) /= n .or. size(problem%upper) /= n) then
         message = 'parameter bounds must have size(problem%start)'
         return
      end if
      if (any(.not. ieee_is_finite(problem%lower)) .or. &
          any(.not. ieee_is_finite(problem%upper))) then
         where (.not. ieee_is_finite(problem%lower)) problem%lower = -default_bound
         where (.not. ieee_is_finite(problem%upper)) problem%upper = default_bound
      end if
      if (any(problem%lower > problem%upper)) then
         message = 'lower parameter bounds exceed upper bounds'
         return
      end if
      problem%start = max(problem%lower, min(problem%upper, problem%start))

      if (problem%raw_n_eq == 0) problem%raw_n_eq = problem%n_eq
      if (problem%raw_n_ineq == 0) problem%raw_n_ineq = problem%n_ineq
      if (problem%n_eq < 0 .or. problem%n_ineq < 0) then
         message = 'constraint dimensions cannot be negative'
         return
      end if
      if (problem%n_eq > 0 .and. .not. associated(problem%eq_fn)) then
         message = 'positive n_eq requires an equality callback'
         return
      end if
      if (problem%n_ineq > 0 .and. .not. associated(problem%ineq_fn)) then
         message = 'positive n_ineq requires an inequality callback'
         return
      end if

      if (problem%n_eq > 0) then
         if (.not. allocated(problem%eq_b)) then
            allocate(problem%eq_b(problem%n_eq))
            problem%eq_b = 0.0_dp
         end if
         if (size(problem%eq_b) /= problem%n_eq) then
            message = 'eq_b must have n_eq elements'
            return
         end if
      else
         if (.not. allocated(problem%eq_b)) allocate(problem%eq_b(0))
      end if

      if (problem%n_ineq > 0) then
         if (.not. allocated(problem%ineq_lower)) then
            allocate(problem%ineq_lower(problem%n_ineq))
            problem%ineq_lower = -default_bound
         end if
         if (.not. allocated(problem%ineq_upper)) then
            allocate(problem%ineq_upper(problem%n_ineq))
            problem%ineq_upper = default_bound
         end if
         if (size(problem%ineq_lower) /= problem%n_ineq .or. &
             size(problem%ineq_upper) /= problem%n_ineq) then
            message = 'inequality bounds must have n_ineq elements'
            return
         end if
         where (.not. ieee_is_finite(problem%ineq_lower)) problem%ineq_lower = -default_bound
         where (.not. ieee_is_finite(problem%ineq_upper)) problem%ineq_upper = default_bound
         if (any(problem%ineq_lower > problem%ineq_upper)) then
            message = 'lower inequality bounds exceed upper bounds'
            return
         end if
      else
         if (.not. allocated(problem%ineq_lower)) allocate(problem%ineq_lower(0))
         if (.not. allocated(problem%ineq_upper)) allocate(problem%ineq_upper(0))
      end if

      status = solnp_success
      message = 'success'
   end subroutine prepare_problem

   subroutine solnp_standardize_problem(problem, standardized, status, message)
      type(solnp_problem), intent(in) :: problem
      type(solnp_problem), intent(out) :: standardized
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message

      integer :: i, m, stat
      character(len=160) :: msg

      standardized = problem
      call prepare_problem(standardized, stat, msg)
      if (stat /= solnp_success) then
         if (present(status)) status = stat
         if (present(message)) message = trim(msg)
         return
      end if
      standardized%standard_form = .true.
      standardized%raw_n_eq = problem%n_eq
      standardized%raw_n_ineq = problem%n_ineq

      if (allocated(standardized%standard_eq_shift)) deallocate(standardized%standard_eq_shift)
      allocate(standardized%standard_eq_shift(problem%n_eq))
      if (problem%n_eq > 0) standardized%standard_eq_shift = standardized%eq_b
      standardized%eq_b = 0.0_dp

      if (allocated(standardized%standard_ineq_lower)) deallocate(standardized%standard_ineq_lower)
      if (allocated(standardized%standard_ineq_upper)) deallocate(standardized%standard_ineq_upper)
      allocate(standardized%standard_ineq_lower(problem%n_ineq))
      allocate(standardized%standard_ineq_upper(problem%n_ineq))
      standardized%standard_ineq_lower = problem%ineq_lower
      standardized%standard_ineq_upper = problem%ineq_upper
      m = 0
      do i = 1, problem%n_ineq
         if (problem%ineq_lower(i) > -0.5_dp * default_bound) m = m + 1
         if (problem%ineq_upper(i) < 0.5_dp * default_bound) m = m + 1
      end do
      standardized%n_ineq = m
      if (allocated(standardized%ineq_lower)) deallocate(standardized%ineq_lower)
      if (allocated(standardized%ineq_upper)) deallocate(standardized%ineq_upper)
      allocate(standardized%ineq_lower(m), standardized%ineq_upper(m))
      standardized%ineq_lower = -default_bound
      standardized%ineq_upper = 0.0_dp

      if (present(status)) status = solnp_success
      if (present(message)) message = 'success'
   end subroutine solnp_standardize_problem

end module rsolnp_problem
