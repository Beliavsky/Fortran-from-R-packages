! SPDX-License-Identifier: LGPL-3.0-only
module pso_benchmarks
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use pso_kinds, only : dp
   use pso_types, only : pso_control, pso_result, pso_objective, pso_gradient, pso_unlimited
   use pso_core, only : psoptim
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   type, public :: pso_test_problem
      character(len=:), allocatable :: name
      procedure(pso_objective), pointer, nopass :: f => null()
      procedure(pso_gradient), pointer, nopass :: grad => null()
      integer :: n = 0
      integer :: maxf = 0
      real(dp) :: objective = 0.0_dp
      integer :: ntest = 100
      real(dp), allocatable :: lower(:), upper(:)
   end type pso_test_problem

   type, public :: pso_test_summary
      type(pso_test_problem) :: problem
      type(pso_result), allocatable :: results(:)
      real(dp) :: cpu_seconds = 0.0_dp
      real(dp) :: elapsed_seconds = 0.0_dp
   end type pso_test_summary

   type, public :: pso_test_statistics
      real(dp) :: objective_mean = 0.0_dp
      real(dp) :: objective_sd = 0.0_dp
      real(dp) :: objective_min = 0.0_dp
      real(dp) :: objective_max = 0.0_dp
      real(dp) :: success_rate = 0.0_dp
      real(dp) :: efficiency = 0.0_dp
      real(dp) :: cpu_seconds = 0.0_dp
      real(dp) :: elapsed_seconds = 0.0_dp
   end type pso_test_statistics

   public :: make_test_problem, run_test_problem, get_success_curve
   public :: test_efficiency, summarize_test
   public :: parabola, parabola_grad, griewank, griewank_grad
   public :: rosenbrock_shifted, rosenbrock_shifted_grad
   public :: rastrigin, rastrigin_grad, ackley, ackley_grad

contains

   function parabola(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = sum(x * x)
   end function parabola

   subroutine parabola_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g = 2.0_dp * x
   end subroutine parabola_grad

   function griewank(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f, p
      integer :: i

      p = 1.0_dp
      do i = 1, size(x)
         p = p * cos(x(i) / sqrt(real(i, dp)))
      end do
      f = sum(x * x) / 4000.0_dp - p + 1.0_dp
   end function griewank

   subroutine griewank_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      real(dp), allocatable :: z(:), dz(:)
      real(dp) :: p
      integer :: i, j

      allocate(z(size(x)), dz(size(x)))
      do i = 1, size(x)
         z(i) = cos(x(i) / sqrt(real(i, dp)))
         dz(i) = -sin(x(i) / sqrt(real(i, dp))) / sqrt(real(i, dp))
      end do
      do i = 1, size(x)
         p = 1.0_dp
         do j = 1, size(x)
            if (j /= i) p = p * z(j)
         end do
         g(i) = x(i) / 2000.0_dp - dz(i) * p
      end do
   end subroutine griewank_grad

   function rosenbrock_shifted(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      real(dp), allocatable :: t0(:), t1(:)
      integer :: n

      n = size(x)
      if (n < 2) then
         f = 0.0_dp
         return
      end if
      allocate(t0(n), t1(n-1))
      t0 = x + 1.0_dp
      t1 = t0(2:n) - t0(1:n-1) * t0(1:n-1)
      f = 100.0_dp * sum(t1 * t1) + sum(x(1:n-1) * x(1:n-1))
   end function rosenbrock_shifted

   subroutine rosenbrock_shifted_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      real(dp), allocatable :: t0(:), t1(:)
      integer :: n

      n = size(x)
      g = 0.0_dp
      if (n < 2) return
      allocate(t0(n), t1(n-1))
      t0 = x + 1.0_dp
      t1 = t0(2:n) - t0(1:n-1) * t0(1:n-1)
      ! This follows the upstream R gradient literally, including 2*t0.
      g(1:n-1) = -400.0_dp * t1 * t0(1:n-1) + 2.0_dp * t0(1:n-1)
      g(2:n) = g(2:n) + 200.0_dp * t1
   end subroutine rosenbrock_shifted_grad

   function rastrigin(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = 10.0_dp * real(size(x), dp) + sum(x*x - 10.0_dp*cos(2.0_dp*pi*x))
   end function rastrigin

   subroutine rastrigin_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g = 2.0_dp*x + 20.0_dp*pi*sin(2.0_dp*pi*x)
   end subroutine rastrigin_grad

   function ackley(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f, n
      n = real(size(x), dp)
      f = -20.0_dp * exp(-0.2_dp * sqrt(sum(x*x)/n)) &
         - exp(sum(cos(2.0_dp*pi*x))/n) + 20.0_dp + exp(1.0_dp)
   end function ackley

   subroutine ackley_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      real(dp) :: n, ss, first, second

      n = real(size(x), dp)
      ss = sum(x*x)
      if (ss > 0.0_dp) then
         first = 4.0_dp * exp(-0.2_dp * sqrt(ss/n)) / sqrt(ss*n)
      else
         first = 0.0_dp
      end if
      second = exp(sum(cos(2.0_dp*pi*x))/n) * 2.0_dp*pi/n
      g = first*x + second*sin(2.0_dp*pi*x)
   end subroutine ackley_grad

   subroutine make_test_problem(name, problem, ntest, dim, maxf, objective, lower, upper)
      character(len=*), intent(in) :: name
      type(pso_test_problem), intent(out) :: problem
      integer, intent(in), optional :: ntest, dim, maxf
      real(dp), intent(in), optional :: objective
      real(dp), intent(in), optional :: lower(:), upper(:)
      integer :: n
      real(dp) :: lo, hi

      problem%name = trim(name)
      select case (trim(name))
      case ("parabola")
         n = 30; problem%maxf = 50000; problem%objective = 1.0e-8_dp
         lo = -100.0_dp; hi = 100.0_dp
         problem%f => parabola; problem%grad => parabola_grad
      case ("griewank")
         n = 2; problem%maxf = 30000; problem%objective = 1.0e-3_dp
         lo = -600.0_dp; hi = 600.0_dp
         problem%f => griewank; problem%grad => griewank_grad
      case ("rosenbrock")
         n = 20; problem%maxf = 100000; problem%objective = 1.0e-4_dp
         lo = -10.0_dp; hi = 10.0_dp
         problem%f => rosenbrock_shifted; problem%grad => rosenbrock_shifted_grad
      case ("rastrigin")
         n = 2; problem%maxf = 3000; problem%objective = 0.0_dp
         lo = -5.12_dp; hi = 5.12_dp
         problem%f => rastrigin; problem%grad => rastrigin_grad
      case ("ackley")
         n = 10; problem%maxf = 5000; problem%objective = 1.0e-4_dp
         lo = -32.0_dp; hi = 32.0_dp
         problem%f => ackley; problem%grad => ackley_grad
      case default
         error stop "make_test_problem: unknown problem name"
      end select

      if (present(dim)) n = dim
      problem%n = n
      problem%ntest = 100
      if (present(ntest)) problem%ntest = ntest
      if (present(maxf)) problem%maxf = maxf
      if (present(objective)) problem%objective = objective
      allocate(problem%lower(n), problem%upper(n))
      problem%lower = lo
      problem%upper = hi
      if (present(lower)) then
         if (size(lower) /= n) error stop "make_test_problem: lower has wrong size"
         problem%lower = lower
      end if
      if (present(upper)) then
         if (size(upper) /= n) error stop "make_test_problem: upper has wrong size"
         problem%upper = upper
      end if
   end subroutine make_test_problem

   subroutine run_test_problem(problem, summary, control)
      type(pso_test_problem), intent(in) :: problem
      type(pso_test_summary), intent(out) :: summary
      type(pso_control), intent(in), optional :: control
      type(pso_control) :: con
      real(dp), allocatable :: par(:)
      real(dp) :: t0, t1
      integer :: i
      integer :: clock0, clock1, clock_rate

      if (.not. associated(problem%f)) error stop "run_test_problem: objective is not associated"
      con = pso_control()
      if (present(control)) con = control
      con%maxf = problem%maxf
      if (con%maxit == 1000) con%maxit = pso_unlimited
      con%abstol = problem%objective

      summary%problem = problem
      allocate(summary%results(problem%ntest), par(problem%n))
      par = ieee_value(0.0_dp, ieee_quiet_nan)
      call cpu_time(t0)
      call system_clock(clock0, clock_rate)
      do i = 1, problem%ntest
         if (associated(problem%grad)) then
            call psoptim(par, problem%f, problem%lower, problem%upper, &
               summary%results(i), con, problem%grad)
         else
            call psoptim(par, problem%f, problem%lower, problem%upper, &
               summary%results(i), con)
         end if
      end do
      call cpu_time(t1)
      call system_clock(clock1)
      summary%cpu_seconds = t1 - t0
      if (clock_rate > 0) then
         summary%elapsed_seconds = real(clock1 - clock0, dp) / real(clock_rate, dp)
      end if
   end subroutine run_test_problem

   function summarize_test(summary) result(stats)
      type(pso_test_summary), intent(in) :: summary
      type(pso_test_statistics) :: stats
      real(dp), allocatable :: values(:)
      real(dp) :: mean_value
      integer :: i, n, nsuccess

      n = size(summary%results)
      if (n < 1) return
      allocate(values(n))
      nsuccess = 0
      do i = 1, n
         values(i) = summary%results(i)%value
         if (values(i) <= summary%problem%objective) nsuccess = nsuccess + 1
      end do
      mean_value = sum(values) / real(n, dp)
      stats%objective_mean = mean_value
      if (n > 1) stats%objective_sd = sqrt(sum((values - mean_value)**2) / real(n - 1, dp))
      stats%objective_min = minval(values)
      stats%objective_max = maxval(values)
      stats%success_rate = real(nsuccess, dp) / real(n, dp)
      stats%efficiency = test_efficiency(summary)
      stats%cpu_seconds = summary%cpu_seconds
      stats%elapsed_seconds = summary%elapsed_seconds
   end function summarize_test

   subroutine get_success_curve(summary, feval, rate)
      type(pso_test_summary), intent(in) :: summary
      integer, allocatable, intent(out) :: feval(:)
      real(dp), allocatable, intent(out) :: rate(:)
      integer, allocatable :: idx(:), fsort(:), ftmp(:)
      real(dp), allocatable :: rtmp(:)
      logical, allocatable :: success(:)
      integer :: n, i, j, ns, nkeep
      real(dp) :: current_rate

      n = size(summary%results)
      allocate(idx(n), fsort(n), success(n), ftmp(n), rtmp(n))
      do i = 1, n
         idx(i) = i
      end do
      call sort_by_feval(summary, idx)
      ns = 0
      nkeep = 0
      do i = 1, n
         j = idx(i)
         fsort(i) = summary%results(j)%function_evaluations
         success(i) = summary%results(j)%value <= summary%problem%objective
         if (success(i)) ns = ns + 1
         current_rate = real(ns, dp) / real(n, dp)
         if (i == 1) then
            nkeep = nkeep + 1
            ftmp(nkeep) = fsort(i)
            rtmp(nkeep) = current_rate
         else if (fsort(i) /= fsort(i-1) .or. &
                  abs(current_rate - rtmp(nkeep)) > tiny(1.0_dp)) then
            nkeep = nkeep + 1
            ftmp(nkeep) = fsort(i)
            rtmp(nkeep) = current_rate
         end if
      end do
      allocate(feval(nkeep), rate(nkeep))
      feval = ftmp(1:nkeep)
      rate = rtmp(1:nkeep)
   end subroutine get_success_curve

   real(dp) function test_efficiency(summary) result(eff)
      type(pso_test_summary), intent(in) :: summary
      integer, allocatable :: feval(:)
      real(dp), allocatable :: rate(:)
      integer :: i, n
      real(dp) :: area

      call get_success_curve(summary, feval, rate)
      n = size(feval)
      area = 0.0_dp
      if (n > 0) then
         area = 0.5_dp * rate(1) * real(feval(1), dp)
         do i = 2, n
            area = area + 0.5_dp * (rate(i) + rate(i-1)) * &
               real(feval(i) - feval(i-1), dp)
         end do
         area = area + 0.5_dp * (1.0_dp + rate(n)) * &
            real(summary%problem%maxf - feval(n), dp)
      else
         area = 0.5_dp * real(summary%problem%maxf, dp)
      end if
      eff = area / real(summary%problem%maxf, dp)
   end function test_efficiency

   subroutine sort_by_feval(summary, idx)
      type(pso_test_summary), intent(in) :: summary
      integer, intent(inout) :: idx(:)
      integer :: i, j, key

      do i = 2, size(idx)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (summary%results(idx(j))%function_evaluations <= &
                summary%results(key)%function_evaluations) exit
            idx(j+1) = idx(j)
            j = j - 1
         end do
         idx(j+1) = key
      end do
   end subroutine sort_by_feval

end module pso_benchmarks
