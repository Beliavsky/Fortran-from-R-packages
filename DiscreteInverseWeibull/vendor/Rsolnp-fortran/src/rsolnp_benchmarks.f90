! SPDX-License-Identifier: GPL-2.0-only
module rsolnp_benchmarks
   use rsolnp_kinds, only : dp
   use rsolnp_types, only : solnp_problem, problem_table_entry, solnp_success, solnp_invalid_problem
   implicit none
   private

   integer, parameter :: id_hs01 = 1, id_hs05 = 5, id_hs06 = 6, id_hs11 = 11
   integer, parameter :: id_hs14 = 14, id_hs24 = 24, id_hs38 = 38, id_hs64 = 64
   integer, parameter :: id_hs65 = 65, id_alkylation = 101, id_box = 102
   integer, parameter :: id_entropy = 103, id_himmelblau5 = 104, id_powell = 105
   integer, parameter :: id_rosen_suzuki = 106, id_wright4 = 107, id_wright9 = 108, id_garch = 109

   type :: benchmark_context
      integer :: id = 0
      real(dp), allocatable :: series(:)
   end type benchmark_context

   public :: solnp_problem_suite, solnp_problems_table

contains

   subroutine solnp_problem_suite(suite, number, problem, status, message)
      character(len=*), intent(in) :: suite
      integer, intent(in) :: number
      type(solnp_problem), intent(out) :: problem
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message

      character(len=32) :: name
      character(len=:), allocatable :: suite_lc
      integer :: stat
      character(len=160) :: msg

      suite_lc = lower(trim(suite))
      name = ''
      if (suite_lc == 'hock-schittkowski' .or. suite_lc == 'hs') then
         write(name, '(a,i0)') 'hs', number
      else if (suite_lc == 'other') then
         select case (number)
         case (1)
            name = 'alkylation'
         case (2)
            name = 'box'
         case (3)
            name = 'entropy'
         case (4)
            name = 'garch'
         case (5)
            name = 'himmelblau5'
         case (6)
            name = 'powell'
         case (7)
            name = 'rosen_suzuki'
         case (8)
            name = 'wright4'
         case (9)
            name = 'wright9'
         case default
            name = ''
         end select
      end if
      if (len_trim(name) == 0) then
         stat = solnp_invalid_problem
         msg = 'unknown suite or problem number'
      else
         call make_benchmark(trim(name), problem, stat, msg)
      end if
      if (present(status)) status = stat
      if (present(message)) message = trim(msg)
   end subroutine solnp_problem_suite

   subroutine solnp_problems_table(entries)
      type(problem_table_entry), allocatable, intent(out) :: entries(:)
      character(len=32), parameter :: other_names(9) = [character(len=32) :: &
         'alkylation', 'box', 'entropy', 'garch', 'himmelblau5', 'powell', &
         'rosen_suzuki', 'wright4', 'wright9']
      integer, parameter :: hs_numbers(68) = [ &
         1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, &
         13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, &
         25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, &
         37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, &
         49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
         61, 62, 63, 64, 65, 110, 118, 119 ]
      integer :: i, k

      allocate(entries(77))
      k = 0
      do i = 1, size(hs_numbers)
         k = k + 1
         entries(k)%suite = 'Hock-Schittkowski'
         write(entries(k)%name, '(a,i0)') 'hs', hs_numbers(i)
         entries(k)%number = hs_numbers(i)
         entries(k)%implemented = is_implemented(trim(entries(k)%name))
      end do
      do i = 1, size(other_names)
         k = k + 1
         entries(k)%suite = 'Other'
         entries(k)%name = other_names(i)
         entries(k)%number = i
         entries(k)%implemented = is_implemented(trim(entries(k)%name))
      end do
   end subroutine solnp_problems_table

   logical function is_implemented(name) result(value)
      character(len=*), intent(in) :: name
      select case (lower(trim(name)))
      case ('hs01', 'hs1', 'hs05', 'hs5', 'hs06', 'hs6', 'hs11', 'hs14', 'hs24', &
            'hs38', 'hs64', 'hs65', 'alkylation', 'box', 'entropy', 'himmelblau5', &
            'powell', 'rosen_suzuki', 'wright4', 'wright9', 'garch')
         value = .true.
      case default
         value = .false.
      end select
   end function is_implemented

   subroutine make_benchmark(name, problem, status, message)
      character(len=*), intent(in) :: name
      type(solnp_problem), intent(out) :: problem
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      type(benchmark_context) :: context
      character(len=:), allocatable :: key
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp) :: aex, aw, qaw, w7

      key = lower(trim(name))
      problem%name = key
      problem%fn => benchmark_objective
      context%id = 0
      select case (key)
      case ('hs01', 'hs1')
         context%id = id_hs01
         call allocate_problem(problem, 2, 0, 0)
         problem%start = [-2.0_dp, 1.0_dp]
         problem%lower = [-1000.0_dp, -1.5_dp]
         problem%upper = [1000.0_dp, 1000.0_dp]
         problem%best_fn = 0.0_dp
         problem%best_par = [1.0_dp, 1.0_dp]
         problem%gr => benchmark_gradient
      case ('hs05', 'hs5')
         context%id = id_hs05
         call allocate_problem(problem, 2, 0, 0)
         problem%start = 0.0_dp
         problem%lower = [-1.5_dp, -3.0_dp]
         problem%upper = [4.0_dp, 3.0_dp]
         problem%best_fn = -(sqrt(3.0_dp) / 2.0_dp + pi / 3.0_dp)
         problem%best_par = [0.5_dp - pi / 3.0_dp, -0.5_dp - pi / 3.0_dp]
         problem%gr => benchmark_gradient
      case ('hs06', 'hs6')
         context%id = id_hs06
         call allocate_problem(problem, 2, 1, 0)
         problem%start = [-1.2_dp, 1.0_dp]
         problem%lower = -1000.0_dp
         problem%upper = 1000.0_dp
         problem%eq_b = 0.0_dp
         problem%best_fn = 0.0_dp
         problem%best_par = [1.0_dp, 1.0_dp]
         problem%gr => benchmark_gradient
         problem%eq_fn => benchmark_equalities
         problem%eq_jac => benchmark_eq_jacobian
      case ('hs11')
         context%id = id_hs11
         call allocate_problem(problem, 2, 0, 1)
         problem%start = [4.9_dp, 0.1_dp]
         problem%lower = -1000.0_dp
         problem%upper = 1000.0_dp
         problem%ineq_lower = 0.0_dp
         problem%ineq_upper = 1.0e10_dp
         aex = 7.5_dp * sqrt(6.0_dp)
         aw = (sqrt(aex * aex + 1.0_dp) + aex) ** (1.0_dp / 3.0_dp)
         qaw = aw * aw
         problem%best_par = [(aw - 1.0_dp / aw) / sqrt(6.0_dp), &
            (qaw - 2.0_dp + 1.0_dp / qaw) / 6.0_dp]
         problem%best_fn = (problem%best_par(1) - 5.0_dp) ** 2 + &
            problem%best_par(2) ** 2 - 25.0_dp
         problem%gr => benchmark_gradient
         problem%ineq_fn => benchmark_inequalities
         problem%ineq_jac => benchmark_ineq_jacobian
      case ('hs14')
         context%id = id_hs14
         call allocate_problem(problem, 2, 1, 1)
         problem%start = [2.0_dp, 2.0_dp]
         problem%lower = -1000.0_dp
         problem%upper = 1000.0_dp
         problem%eq_b = 0.0_dp
         problem%ineq_lower = 0.0_dp
         problem%ineq_upper = 1.0e8_dp
         w7 = sqrt(7.0_dp)
         problem%best_par = [(w7 - 1.0_dp) / 2.0_dp, (w7 + 1.0_dp) / 4.0_dp]
         problem%best_fn = 9.0_dp - 23.0_dp * w7 / 8.0_dp
         problem%gr => benchmark_gradient
         problem%eq_fn => benchmark_equalities
         problem%eq_jac => benchmark_eq_jacobian
         problem%ineq_fn => benchmark_inequalities
         problem%ineq_jac => benchmark_ineq_jacobian
      case ('hs24')
         context%id = id_hs24
         call allocate_problem(problem, 2, 0, 3)
         problem%start = [500.0_dp, 500.0_dp]
         problem%lower = 0.0_dp
         problem%upper = 1000.0_dp
         problem%ineq_lower = 0.0_dp
         problem%ineq_upper = 1.0e8_dp
         problem%best_fn = -1.0_dp
         problem%best_par = [3.0_dp, sqrt(3.0_dp)]
         problem%gr => benchmark_gradient
         problem%ineq_fn => benchmark_inequalities
         problem%ineq_jac => benchmark_ineq_jacobian
      case ('hs38')
         context%id = id_hs38
         call allocate_problem(problem, 4, 0, 0)
         problem%start = [-3.0_dp, -1.0_dp, -3.0_dp, -1.0_dp]
         problem%lower = -10.0_dp
         problem%upper = 10.0_dp
         problem%best_fn = 0.0_dp
         problem%best_par = 1.0_dp
         problem%gr => benchmark_gradient
      case ('hs64')
         context%id = id_hs64
         call allocate_problem(problem, 3, 0, 1)
         problem%start = 1.0_dp
         problem%lower = 1.0e-5_dp
         problem%upper = 1000.0_dp
         problem%ineq_lower = 0.0_dp
         problem%ineq_upper = 1.0e8_dp
         problem%best_fn = 6299.84242821_dp
         problem%best_par = [108.734717597_dp, 85.1261394257_dp, 204.324707858_dp]
         problem%gr => benchmark_gradient
         problem%ineq_fn => benchmark_inequalities
         problem%ineq_jac => benchmark_ineq_jacobian
      case ('hs65')
         context%id = id_hs65
         call allocate_problem(problem, 3, 0, 1)
         problem%start = [-5.0_dp, 5.0_dp, 0.0_dp]
         problem%lower = [-4.5_dp, -4.5_dp, -5.0_dp]
         problem%upper = [4.5_dp, 4.5_dp, 5.0_dp]
         problem%ineq_lower = 0.0_dp
         problem%ineq_upper = 1.0e8_dp
         problem%best_fn = 0.953528856757_dp
         problem%best_par = [3.65046182158_dp, 3.65046168940_dp, 4.62041750754_dp]
         problem%gr => benchmark_gradient
         problem%ineq_fn => benchmark_inequalities
         problem%ineq_jac => benchmark_ineq_jacobian
      case ('alkylation')
         context%id = id_alkylation
         call allocate_problem(problem, 10, 3, 4)
         problem%start = [17.45_dp, 12.0_dp, 110.0_dp, 30.0_dp, 19.74_dp, &
            89.2_dp, 92.8_dp, 8.0_dp, 3.6_dp, 155.0_dp]
         problem%lower = [0.0_dp, 0.0_dp, 0.0_dp, 10.0_dp, 0.0_dp, 85.0_dp, &
            10.0_dp, 3.0_dp, 1.0_dp, 145.0_dp]
         problem%upper = [20.0_dp, 16.0_dp, 120.0_dp, 50.0_dp, 20.0_dp, 93.0_dp, &
            95.0_dp, 12.0_dp, 4.0_dp, 162.0_dp]
         problem%eq_b = 0.0_dp
         problem%ineq_lower = [0.99_dp, 0.99_dp, 0.9_dp, 0.99_dp]
         problem%ineq_upper = [100.0_dp / 99.0_dp, 100.0_dp / 99.0_dp, &
            10.0_dp / 9.0_dp, 100.0_dp / 99.0_dp]
         problem%best_fn = -172.642_dp
         problem%best_par = [16.996427_dp, 16.0_dp, 57.685751_dp, 30.324940_dp, &
            20.0_dp, 90.565147_dp, 95.0_dp, 10.590461_dp, 1.561636_dp, 153.535354_dp]
         problem%eq_fn => benchmark_equalities
         problem%ineq_fn => benchmark_inequalities
      case ('box')
         context%id = id_box
         call allocate_problem(problem, 3, 1, 0)
         problem%start = [1.1_dp, 1.1_dp, 9.0_dp]
         problem%lower = 1.0_dp
         problem%upper = 10.0_dp
         problem%eq_b = 100.0_dp
         problem%best_fn = -48.11252_dp
         problem%best_par = [2.886751_dp, 2.886751_dp, 5.773503_dp]
         problem%gr => benchmark_gradient
         problem%eq_fn => benchmark_equalities
         problem%eq_jac => benchmark_eq_jacobian
      case ('entropy')
         context%id = id_entropy
         call allocate_problem(problem, 10, 1, 0)
         problem%start = 1.0_dp
         problem%lower = 1.0e-8_dp
         problem%upper = 1000.0_dp
         problem%eq_b = 10.0_dp
         problem%best_fn = 0.1854782_dp
         problem%best_par = [2.2801555_dp, 0.8577605_dp, 0.8577605_dp, &
            0.8577605_dp, 0.8577605_dp, 0.8577605_dp, 0.8577605_dp, &
            0.8577605_dp, 0.8577605_dp, 0.8577605_dp]
         problem%gr => benchmark_gradient
         problem%eq_fn => benchmark_equalities
         problem%eq_jac => benchmark_eq_jacobian
      case ('garch')
         context%id = id_garch
         call allocate_problem(problem, 4, 0, 1)
         problem%start = [0.0_dp, 0.1_dp, 0.05_dp, 0.9_dp]
         problem%lower = [-1.0_dp, 1.0e-12_dp, 1.0e-12_dp, 1.0e-12_dp]
         problem%upper = [1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp]
         problem%ineq_lower = -1.0_dp
         problem%ineq_upper = 0.0_dp
         problem%best_fn = 1074.36_dp
         problem%best_par = [-0.006184353_dp, 0.010760430_dp, 0.153408326_dp, 0.805877422_dp]
         allocate(context%series(1500))
         call generate_garch_series(context%series)
         problem%gr => benchmark_gradient
         problem%ineq_fn => benchmark_inequalities
         problem%ineq_jac => benchmark_ineq_jacobian
      case ('himmelblau5')
         context%id = id_himmelblau5
         call allocate_problem(problem, 2, 1, 1)
         problem%start = [1.0_dp, 1.0_dp]
         problem%lower = -5.0_dp
         problem%upper = 5.0_dp
         problem%eq_b = 0.0_dp
         problem%ineq_lower = -1000.0_dp
         problem%ineq_upper = 25.0_dp
         problem%best_fn = 11.7544_dp
         problem%best_par = [4.74342_dp, 1.58114_dp]
         problem%gr => benchmark_gradient
         problem%eq_fn => benchmark_equalities
         problem%eq_jac => benchmark_eq_jacobian
         problem%ineq_fn => benchmark_inequalities
         problem%ineq_jac => benchmark_ineq_jacobian
      case ('powell')
         context%id = id_powell
         call allocate_problem(problem, 5, 3, 0)
         problem%start = [-2.0_dp, 2.0_dp, 2.0_dp, -1.0_dp, -1.0_dp]
         problem%lower = -10.0_dp
         problem%upper = 10.0_dp
         problem%eq_b = [10.0_dp, 0.0_dp, -1.0_dp]
         problem%best_fn = 0.05394985_dp
         problem%best_par = [-1.717144_dp, 1.595710_dp, 1.827245_dp, 0.763643_dp, 0.763643_dp]
         problem%gr => benchmark_gradient
         problem%eq_fn => benchmark_equalities
         problem%eq_jac => benchmark_eq_jacobian
      case ('rosen_suzuki')
         context%id = id_rosen_suzuki
         call allocate_problem(problem, 4, 0, 3)
         problem%start = 1.0_dp
         problem%lower = -10.0_dp
         problem%upper = 10.0_dp
         problem%ineq_lower = 0.0_dp
         problem%ineq_upper = 1000.0_dp
         problem%best_fn = -44.0_dp
         problem%best_par = [2.502771e-7_dp, 0.9999997_dp, 2.0_dp, -1.0_dp]
         problem%gr => benchmark_gradient
         problem%ineq_fn => benchmark_inequalities
         problem%ineq_jac => benchmark_ineq_jacobian
      case ('wright4')
         context%id = id_wright4
         call allocate_problem(problem, 5, 3, 0)
         problem%start = 1.0_dp
         problem%lower = -10.0_dp
         problem%upper = 10.0_dp
         problem%eq_b = [2.0_dp + 3.0_dp * sqrt(2.0_dp), &
            -2.0_dp + 2.0_dp * sqrt(2.0_dp), 2.0_dp]
         problem%best_fn = 0.02931083_dp
         problem%best_par = [1.116635_dp, 1.220442_dp, 1.537785_dp, 1.972769_dp, 1.791096_dp]
         problem%gr => benchmark_gradient
         problem%eq_fn => benchmark_equalities
         problem%eq_jac => benchmark_eq_jacobian
      case ('wright9')
         context%id = id_wright9
         call allocate_problem(problem, 5, 0, 3)
         problem%start = 1.0_dp
         problem%lower = -5.0_dp
         problem%upper = 5.0_dp
         problem%ineq_lower = [-100.0_dp, -2.0_dp, 5.0_dp]
         problem%ineq_upper = [20.0_dp, 100.0_dp, 100.0_dp]
         problem%best_fn = -210.4078_dp
         problem%best_par = [-0.08145219_dp, 3.69237756_dp, 2.48741102_dp, &
            0.37713392_dp, 0.17398257_dp]
         problem%gr => benchmark_gradient
         problem%ineq_fn => benchmark_inequalities
         problem%ineq_jac => benchmark_ineq_jacobian
      case default
         status = solnp_invalid_problem
         message = 'benchmark definition is not translated'
         return
      end select
      allocate(problem%data, source=context)
      status = solnp_success
      message = 'success'
   end subroutine make_benchmark

   subroutine allocate_problem(problem, n, n_eq, n_ineq)
      type(solnp_problem), intent(inout) :: problem
      integer, intent(in) :: n, n_eq, n_ineq

      problem%n = n
      problem%n_eq = n_eq
      problem%n_ineq = n_ineq
      problem%raw_n_eq = n_eq
      problem%raw_n_ineq = n_ineq
      allocate(problem%start(n), problem%lower(n), problem%upper(n), problem%best_par(n))
      allocate(problem%eq_b(n_eq), problem%ineq_lower(n_ineq), problem%ineq_upper(n_ineq))
      problem%start = 0.0_dp
      problem%lower = -1.0e20_dp
      problem%upper = 1.0e20_dp
      problem%eq_b = 0.0_dp
      problem%ineq_lower = -1.0e20_dp
      problem%ineq_upper = 1.0e20_dp
      problem%best_par = 0.0_dp
   end subroutine allocate_problem

   subroutine benchmark_objective(x, value, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: data
      integer :: id
      real(dp) :: a, norm_diff, exponent_limit

      id = get_id(data)
      select case (id)
      case (id_hs01)
         value = 100.0_dp * (x(2) - x(1) ** 2) ** 2 + (1.0_dp - x(1)) ** 2
      case (id_hs05)
         value = sin(x(1) + x(2)) + (x(1) - x(2)) ** 2 - 1.5_dp * x(1) + 2.5_dp * x(2) + 1.0_dp
      case (id_hs06)
         value = (1.0_dp - x(1)) ** 2
      case (id_hs11)
         value = (x(1) - 5.0_dp) ** 2 + x(2) ** 2 - 25.0_dp
      case (id_hs14)
         value = (x(1) - 2.0_dp) ** 2 + (x(2) - 1.0_dp) ** 2
      case (id_hs24)
         a = sqrt(3.0_dp)
         value = ((x(1) - 3.0_dp) ** 2 - 9.0_dp) * x(2) ** 3 / (27.0_dp * a)
      case (id_hs38)
         value = 100.0_dp * (x(2) - x(1) ** 2) ** 2 + (1.0_dp - x(1)) ** 2 + &
            90.0_dp * (x(4) - x(3) ** 2) ** 2 + (1.0_dp - x(3)) ** 2 + &
            10.1_dp * ((x(2) - 1.0_dp) ** 2 + (x(4) - 1.0_dp) ** 2) + &
            19.8_dp * (x(2) - 1.0_dp) * (x(4) - 1.0_dp)
      case (id_hs64)
         value = 5.0_dp * x(1) + 5.0e4_dp / x(1) + 20.0_dp * x(2) + &
            7.2e4_dp / x(2) + 10.0_dp * x(3) + 1.44e5_dp / x(3)
      case (id_hs65)
         value = (x(1) - x(2)) ** 2 + ((x(1) + x(2) - 10.0_dp) / 3.0_dp) ** 2 + &
            (x(3) - 5.0_dp) ** 2
      case (id_alkylation)
         value = -0.63_dp * x(4) * x(7) + 50.4_dp * x(1) + 3.5_dp * x(2) + x(3) + 33.6_dp * x(5)
      case (id_box)
         value = -x(1) * x(2) * x(3)
      case (id_entropy)
         norm_diff = sqrt(sum((x - 1.0_dp) ** 2))
         value = -sum(log(x)) - log(norm_diff + 0.1_dp)
      case (id_garch)
         call garch_objective(x, data, value)
      case (id_himmelblau5)
         value = (x(1) - 5.0_dp) ** 2 + (x(2) - 5.0_dp) ** 2
      case (id_powell)
         exponent_limit = log(huge(1.0_dp)) - 20.0_dp
         if (product(x) >= exponent_limit) then
            value = exp(exponent_limit)
         else if (product(x) <= log(tiny(1.0_dp))) then
            value = 0.0_dp
         else
            value = exp(product(x))
         end if
      case (id_rosen_suzuki)
         value = x(1) ** 2 + x(2) ** 2 + 2.0_dp * x(3) ** 2 + x(4) ** 2 - &
            5.0_dp * x(1) - 5.0_dp * x(2) - 21.0_dp * x(3) + 7.0_dp * x(4)
      case (id_wright4)
         value = (x(1) - 1.0_dp) ** 2 + (x(1) - x(2)) ** 2 + (x(2) - x(3)) ** 3 + &
            (x(3) - x(4)) ** 4 + (x(4) - x(5)) ** 4
      case (id_wright9)
         value = 10.0_dp * x(1) * x(4) - 6.0_dp * x(3) * x(2) ** 2 + &
            x(2) * x(1) ** 3 + 9.0_dp * sin(x(5) - x(3)) + x(5) ** 4 * x(4) ** 2 * x(2) ** 3
      case default
         value = huge(1.0_dp)
      end select
   end subroutine benchmark_objective

   subroutine benchmark_gradient(x, gradient, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
      class(*), intent(in), optional :: data
      integer :: id, i, j
      real(dp) :: a, fval, prod_term, prod_other, v1, v2, norm_diff

      id = get_id(data)
      gradient = 0.0_dp
      select case (id)
      case (id_hs01)
         gradient(2) = 200.0_dp * (x(2) - x(1) ** 2)
         gradient(1) = -2.0_dp * (x(1) * (gradient(2) - 1.0_dp) + 1.0_dp)
      case (id_hs05)
         gradient(1) = cos(x(1) + x(2)) + 2.0_dp * (x(1) - x(2)) - 1.5_dp
         gradient(2) = cos(x(1) + x(2)) - 2.0_dp * (x(1) - x(2)) + 2.5_dp
      case (id_hs06)
         gradient = [-2.0_dp * (1.0_dp - x(1)), 0.0_dp]
      case (id_hs11)
         gradient = [2.0_dp * (x(1) - 5.0_dp), 2.0_dp * x(2)]
      case (id_hs14)
         gradient = [2.0_dp * (x(1) - 2.0_dp), 2.0_dp * (x(2) - 1.0_dp)]
      case (id_hs24)
         a = sqrt(3.0_dp)
         gradient(1) = 2.0_dp * (x(1) - 3.0_dp) * x(2) ** 3 / (27.0_dp * a)
         gradient(2) = ((x(1) - 3.0_dp) ** 2 - 9.0_dp) * x(2) ** 2 / (9.0_dp * a)
      case (id_hs38)
         gradient(1) = -400.0_dp * x(1) * (x(2) - x(1) ** 2) - 2.0_dp * (1.0_dp - x(1))
         gradient(2) = 200.0_dp * (x(2) - x(1) ** 2) + 20.2_dp * (x(2) - 1.0_dp) + &
            19.8_dp * (x(4) - 1.0_dp)
         gradient(3) = -360.0_dp * x(3) * (x(4) - x(3) ** 2) - 2.0_dp * (1.0_dp - x(3))
         gradient(4) = 180.0_dp * (x(4) - x(3) ** 2) + 20.2_dp * (x(4) - 1.0_dp) + &
            19.8_dp * (x(2) - 1.0_dp)
      case (id_hs64)
         gradient = [5.0_dp - 5.0e4_dp / x(1) ** 2, 20.0_dp - 7.2e4_dp / x(2) ** 2, &
            10.0_dp - 1.44e5_dp / x(3) ** 2]
      case (id_hs65)
         v1 = 2.0_dp * (x(1) - x(2))
         v2 = 2.0_dp * (x(1) + x(2) - 10.0_dp) / 9.0_dp
         gradient = [v1 + v2, -v1 + v2, 2.0_dp * (x(3) - 5.0_dp)]
      case (id_box)
         gradient = [-x(2) * x(3), -x(1) * x(3), -x(1) * x(2)]
      case (id_entropy)
         norm_diff = sqrt(sum((x - 1.0_dp) ** 2))
         gradient = -1.0_dp / x
         if (norm_diff > tiny(1.0_dp)) gradient = gradient - &
            (x - 1.0_dp) / (norm_diff * (norm_diff + 0.1_dp))
      case (id_garch)
         call garch_gradient(x, data, gradient)
      case (id_himmelblau5)
         gradient = [2.0_dp * (x(1) - 5.0_dp), 2.0_dp * (x(2) - 5.0_dp)]
      case (id_powell)
         prod_term = product(x)
         if (prod_term >= log(huge(1.0_dp)) - 20.0_dp) then
            fval = exp(log(huge(1.0_dp)) - 20.0_dp)
         else if (prod_term <= log(tiny(1.0_dp))) then
            fval = 0.0_dp
         else
            fval = exp(prod_term)
         end if
         do i = 1, size(x)
            if (abs(x(i)) > tiny(1.0_dp)) then
               gradient(i) = fval * prod_term / x(i)
            else
               prod_other = 1.0_dp
               do j = 1, size(x)
                  if (j /= i) prod_other = prod_other * x(j)
               end do
               gradient(i) = fval * prod_other
            end if
         end do
      case (id_rosen_suzuki)
         gradient = [2.0_dp * x(1) - 5.0_dp, 2.0_dp * x(2) - 5.0_dp, &
            4.0_dp * x(3) - 21.0_dp, 2.0_dp * x(4) + 7.0_dp]
      case (id_wright4)
         gradient(1) = 2.0_dp * (x(1) - 1.0_dp) + 2.0_dp * (x(1) - x(2))
         gradient(2) = -2.0_dp * (x(1) - x(2)) + 3.0_dp * (x(2) - x(3)) ** 2
         gradient(3) = -3.0_dp * (x(2) - x(3)) ** 2 + 4.0_dp * (x(3) - x(4)) ** 3
         gradient(4) = -4.0_dp * (x(3) - x(4)) ** 3 + 4.0_dp * (x(4) - x(5)) ** 3
         gradient(5) = -4.0_dp * (x(4) - x(5)) ** 3
      case (id_wright9)
         gradient(1) = 10.0_dp * x(4) + 3.0_dp * x(2) * x(1) ** 2
         gradient(2) = -12.0_dp * x(3) * x(2) + x(1) ** 3 + &
            3.0_dp * x(5) ** 4 * x(4) ** 2 * x(2) ** 2
         gradient(3) = -6.0_dp * x(2) ** 2 - 9.0_dp * cos(x(5) - x(3))
         gradient(4) = 10.0_dp * x(1) + 2.0_dp * x(5) ** 4 * x(4) * x(2) ** 3
         gradient(5) = 9.0_dp * cos(x(5) - x(3)) + 4.0_dp * x(5) ** 3 * x(4) ** 2 * x(2) ** 3
      end select
   end subroutine benchmark_gradient

   subroutine benchmark_equalities(x, value, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:)
      class(*), intent(in), optional :: data
      integer :: id

      id = get_id(data)
      value = 0.0_dp
      select case (id)
      case (id_hs06)
         value(1) = 10.0_dp * (x(2) - x(1) ** 2)
      case (id_hs14)
         value(1) = x(1) - 2.0_dp * x(2) + 1.0_dp
      case (id_alkylation)
         value(1) = 98.0_dp * x(3) - 0.1_dp * x(4) * x(6) * x(9) - x(3) * x(6)
         value(2) = 1000.0_dp * x(2) + 100.0_dp * x(5) - 100.0_dp * x(1) * x(8)
         value(3) = 122.0_dp * x(4) - 100.0_dp * x(1) - 100.0_dp * x(5)
      case (id_box)
         value(1) = 4.0_dp * x(1) * x(2) + 2.0_dp * x(2) * x(3) + 2.0_dp * x(3) * x(1)
      case (id_entropy)
         value(1) = sum(x)
      case (id_himmelblau5)
         value(1) = x(1) - 3.0_dp * x(2)
      case (id_powell)
         value(1) = sum(x ** 2)
         value(2) = x(2) * x(3) - 5.0_dp * x(4) * x(5)
         value(3) = x(1) ** 3 + x(2) ** 3
      case (id_wright4)
         value(1) = x(1) + x(2) ** 2 + x(3) ** 3
         value(2) = x(2) - x(3) ** 2 + x(4)
         value(3) = x(1) * x(5)
      end select
   end subroutine benchmark_equalities

   subroutine benchmark_eq_jacobian(x, jacobian, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: jacobian(:, :)
      class(*), intent(in), optional :: data
      integer :: id

      id = get_id(data)
      jacobian = 0.0_dp
      select case (id)
      case (id_hs06)
         jacobian(1, :) = [-20.0_dp * x(1), 10.0_dp]
      case (id_hs14)
         jacobian(1, :) = [1.0_dp, -2.0_dp]
      case (id_box)
         jacobian(1, :) = [4.0_dp * x(2) + 2.0_dp * x(3), &
            4.0_dp * x(1) + 2.0_dp * x(3), 2.0_dp * x(2) + 2.0_dp * x(1)]
      case (id_entropy)
         jacobian(1, :) = 1.0_dp
      case (id_himmelblau5)
         jacobian(1, :) = [1.0_dp, -3.0_dp]
      case (id_powell)
         jacobian(1, :) = 2.0_dp * x
         jacobian(2, :) = [0.0_dp, x(3), x(2), -5.0_dp * x(5), -5.0_dp * x(4)]
         jacobian(3, :) = [3.0_dp * x(1) ** 2, 3.0_dp * x(2) ** 2, 0.0_dp, 0.0_dp, 0.0_dp]
      case (id_wright4)
         jacobian(1, :) = [1.0_dp, 2.0_dp * x(2), 3.0_dp * x(3) ** 2, 0.0_dp, 0.0_dp]
         jacobian(2, :) = [0.0_dp, 1.0_dp, -2.0_dp * x(3), 1.0_dp, 0.0_dp]
         jacobian(3, :) = [x(5), 0.0_dp, 0.0_dp, 0.0_dp, x(1)]
      end select
   end subroutine benchmark_eq_jacobian

   subroutine benchmark_inequalities(x, value, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:)
      class(*), intent(in), optional :: data
      integer :: id
      real(dp) :: a

      id = get_id(data)
      value = 0.0_dp
      select case (id)
      case (id_hs11)
         value(1) = -x(1) ** 2 + x(2)
      case (id_hs14)
         value(1) = 1.0_dp - 0.25_dp * x(1) ** 2 - x(2) ** 2
      case (id_hs24)
         a = sqrt(3.0_dp)
         value = [x(1) / a - x(2), x(1) + a * x(2), 6.0_dp - a * x(2) - x(1)]
      case (id_hs64)
         value(1) = 1.0_dp - 4.0_dp / x(1) - 32.0_dp / x(2) - 120.0_dp / x(3)
      case (id_hs65)
         value(1) = 48.0_dp - x(1) ** 2 - x(2) ** 2 - x(3) ** 2
      case (id_alkylation)
         value(1) = (1.12_dp * x(1) + 0.13167_dp * x(1) * x(8) - &
            0.00667_dp * x(1) * x(8) ** 2) / x(4)
         value(2) = (1.098_dp * x(8) - 0.038_dp * x(8) ** 2 + 0.325_dp * x(6) + 57.25_dp) / x(7)
         value(3) = (-0.222_dp * x(10) + 35.82_dp) / x(9)
         value(4) = (3.0_dp * x(7) - 133.0_dp) / x(10)
      case (id_garch)
         value(1) = x(3) + x(4) - 1.0_dp
      case (id_himmelblau5)
         value(1) = x(1) ** 2 + x(2) ** 2
      case (id_rosen_suzuki)
         value(1) = 8.0_dp - sum(x ** 2) - x(1) + x(2) - x(3) + x(4)
         value(2) = 10.0_dp - x(1) ** 2 - 2.0_dp * x(2) ** 2 - x(3) ** 2 - &
            2.0_dp * x(4) ** 2 + x(1) + x(4)
         value(3) = 5.0_dp - 2.0_dp * x(1) ** 2 - x(2) ** 2 - x(3) ** 2 - &
            2.0_dp * x(1) + x(2) + x(4)
      case (id_wright9)
         value(1) = sum(x ** 2)
         value(2) = x(1) ** 2 * x(3) - x(4) * x(5)
         value(3) = x(2) ** 2 * x(4) + 10.0_dp * x(1) * x(5)
      end select
   end subroutine benchmark_inequalities

   subroutine benchmark_ineq_jacobian(x, jacobian, data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: jacobian(:, :)
      class(*), intent(in), optional :: data
      integer :: id
      real(dp) :: a

      id = get_id(data)
      jacobian = 0.0_dp
      select case (id)
      case (id_hs11)
         jacobian(1, :) = [-2.0_dp * x(1), 1.0_dp]
      case (id_hs14)
         jacobian(1, :) = [-0.5_dp * x(1), -2.0_dp * x(2)]
      case (id_hs24)
         a = sqrt(3.0_dp)
         jacobian(1, :) = [1.0_dp / a, -1.0_dp]
         jacobian(2, :) = [1.0_dp, a]
         jacobian(3, :) = [-1.0_dp, -a]
      case (id_hs64)
         jacobian(1, :) = [4.0_dp / x(1) ** 2, 32.0_dp / x(2) ** 2, 120.0_dp / x(3) ** 2]
      case (id_hs65)
         jacobian(1, :) = -2.0_dp * x
      case (id_garch)
         jacobian(1, :) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
      case (id_himmelblau5)
         jacobian(1, :) = 2.0_dp * x
      case (id_rosen_suzuki)
         jacobian(1, :) = [-2.0_dp * x(1) - 1.0_dp, -2.0_dp * x(2) + 1.0_dp, &
            -2.0_dp * x(3) - 1.0_dp, -2.0_dp * x(4) + 1.0_dp]
         jacobian(2, :) = [-2.0_dp * x(1) + 1.0_dp, -4.0_dp * x(2), &
            -2.0_dp * x(3), -4.0_dp * x(4) + 1.0_dp]
         jacobian(3, :) = [-4.0_dp * x(1) - 2.0_dp, -2.0_dp * x(2) + 1.0_dp, &
            -2.0_dp * x(3), 1.0_dp]
      case (id_wright9)
         jacobian(1, :) = 2.0_dp * x
         jacobian(2, :) = [2.0_dp * x(1) * x(3), 0.0_dp, x(1) ** 2, -x(5), -x(4)]
         jacobian(3, :) = [10.0_dp * x(5), 2.0_dp * x(2) * x(4), 0.0_dp, x(2) ** 2, 10.0_dp * x(1)]
      end select
   end subroutine benchmark_ineq_jacobian

   subroutine generate_garch_series(series)
      real(dp), intent(out) :: series(:)
      real(dp), allocatable :: sigma2(:)
      real(dp) :: mu, omega, alpha, beta, z
      integer :: t, state

      allocate(sigma2(size(series)))
      state = 100
      mu = -0.006184353_dp
      omega = 0.010760430_dp
      alpha = 0.153408326_dp
      beta = 0.805877422_dp
      series(1) = 0.1315172_dp
      sigma2(1) = 0.2211_dp
      do t = 2, size(series)
         sigma2(t) = omega + alpha * (series(t - 1) - mu) ** 2 + beta * sigma2(t - 1)
         call normal_random(state, z)
         series(t) = mu + z * sqrt(sigma2(t))
      end do
   end subroutine generate_garch_series

   subroutine normal_random(state, value)
      integer, intent(inout) :: state
      real(dp), intent(out) :: value
      real(dp) :: u1, u2
      real(dp), parameter :: pi = acos(-1.0_dp)

      call uniform_random(state, u1)
      call uniform_random(state, u2)
      u1 = max(u1, tiny(1.0_dp))
      value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
   end subroutine normal_random

   subroutine uniform_random(state, value)
      integer, intent(inout) :: state
      real(dp), intent(out) :: value
      integer :: k

      k = state / 127773
      state = 16807 * (state - k * 127773) - 2836 * k
      if (state <= 0) state = state + 2147483647
      value = real(state, dp) / 2147483647.0_dp
   end subroutine uniform_random

   subroutine garch_objective(par, data, value)
      real(dp), intent(in) :: par(:)
      class(*), intent(in), optional :: data
      real(dp), intent(out) :: value
      real(dp), allocatable :: residual(:), sigma2(:)
      real(dp) :: mu, omega, alpha, beta
      integer :: t, n

      value = huge(1.0_dp)
      if (.not. present(data)) return
      select type (context => data)
      type is (benchmark_context)
         if (.not. allocated(context%series)) return
         n = size(context%series)
         allocate(residual(n), sigma2(n))
         mu = par(1)
         omega = par(2)
         alpha = par(3)
         beta = par(4)
         residual = context%series - mu
         sigma2(1) = sum(residual ** 2) / real(n, dp)
         do t = 2, n
            sigma2(t) = omega + alpha * residual(t - 1) ** 2 + beta * sigma2(t - 1)
            if (sigma2(t) <= tiny(1.0_dp)) return
         end do
         value = 0.5_dp * sum(log(2.0_dp * acos(-1.0_dp)) + log(sigma2) + residual ** 2 / sigma2)
      class default
         return
      end select
   end subroutine garch_objective

   subroutine garch_gradient(par, data, gradient)
      real(dp), intent(in) :: par(:)
      class(*), intent(in), optional :: data
      real(dp), intent(out) :: gradient(:)
      real(dp), allocatable :: residual(:), sigma2(:), dmu(:), domega(:), dalpha(:), dbeta(:), weight(:)
      real(dp) :: mu, omega, alpha, beta
      integer :: t, n

      gradient = 0.0_dp
      if (.not. present(data)) return
      select type (context => data)
      type is (benchmark_context)
         if (.not. allocated(context%series)) return
         n = size(context%series)
         allocate(residual(n), sigma2(n), dmu(n), domega(n), dalpha(n), dbeta(n), weight(n))
         mu = par(1)
         omega = par(2)
         alpha = par(3)
         beta = par(4)
         residual = context%series - mu
         sigma2(1) = sum(residual ** 2) / real(n, dp)
         dmu = 0.0_dp
         domega = 0.0_dp
         dalpha = 0.0_dp
         dbeta = 0.0_dp
         dmu(1) = -2.0_dp * sum(residual) / real(n, dp)
         do t = 2, n
            domega(t) = 1.0_dp + beta * domega(t - 1)
            dalpha(t) = residual(t - 1) ** 2 + beta * dalpha(t - 1)
            dbeta(t) = sigma2(t - 1) + beta * dbeta(t - 1)
            dmu(t) = -2.0_dp * alpha * residual(t - 1) + beta * dmu(t - 1)
            sigma2(t) = omega + alpha * residual(t - 1) ** 2 + beta * sigma2(t - 1)
            if (sigma2(t) <= tiny(1.0_dp)) then
               gradient = huge(1.0_dp) / 100.0_dp
               return
            end if
         end do
         weight = 0.5_dp * (1.0_dp / sigma2 - residual ** 2 / sigma2 ** 2)
         gradient(1) = -sum(residual / sigma2) + sum(weight * dmu)
         gradient(2) = sum(weight * domega)
         gradient(3) = sum(weight * dalpha)
         gradient(4) = sum(weight * dbeta)
      class default
         return
      end select
   end subroutine garch_gradient

   integer function get_id(data) result(id)
      class(*), intent(in), optional :: data
      id = 0
      if (.not. present(data)) return
      select type (data)
      type is (benchmark_context)
         id = data%id
      class default
         id = 0
      end select
   end function get_id

   pure function lower(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, code
      out = text
      do i = 1, len(text)
         code = iachar(out(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
      end do
   end function lower

end module rsolnp_benchmarks
