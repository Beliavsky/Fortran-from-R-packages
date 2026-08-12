! psoptim-fortran: modern Fortran translation of psoptim 1.0
! Original package copyright: Krzysztof Ciupke
! License: GPL-2.0-or-later. See COPYING and original/DESCRIPTION.
module psoptim
   use, intrinsic :: iso_fortran_env, only : int64
   use psoptim_rng, only : dp, ps_rng
   implicit none
   private

   public :: dp, ps_control, ps_result, ps_optimize
   public :: ps_scalar_objective

   abstract interface
      function ps_scalar_objective(x) result(value)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function ps_scalar_objective
   end interface

   type :: ps_control
      integer :: n = 100
      integer :: max_loop = 100
      real(dp) :: w = 0.9_dp
      real(dp) :: c1 = 0.2_dp
      real(dp) :: c2 = 0.2_dp
      integer(int64) :: seed = 10_int64
      logical :: legacy_initial_bounds = .true.
      logical :: legacy_velocity_initialization = .true.
      logical :: legacy_velocity_clip = .true.
      logical :: store_history = .true.
   end type ps_control

   type :: ps_result
      real(dp), allocatable :: sol(:)
      real(dp) :: val = -huge(1.0_dp)
      integer :: iterations = 0
      integer(int64) :: nfe = 0_int64
      real(dp), allocatable :: best_history(:)
      real(dp), allocatable :: mean_history(:)
      real(dp), allocatable :: final_population(:,:)
      real(dp), allocatable :: final_velocity(:,:)
      real(dp), allocatable :: personal_best(:,:)
   end type ps_result

contains

   subroutine ps_optimize(fun, xmin, xmax, vmax, result, control)
      procedure(ps_scalar_objective) :: fun
      real(dp), intent(in) :: xmin(:), xmax(:), vmax(:)
      type(ps_result), intent(out) :: result
      type(ps_control), intent(in), optional :: control

      type(ps_control) :: ctl
      type(ps_rng) :: rng
      real(dp), allocatable :: x(:,:), v(:,:), pbest(:,:)
      real(dp), allocatable :: fx(:), fpbest(:), gbest(:)
      real(dp) :: gbest_value, candidate_value
      real(dp) :: r1, r2, x_prev, vlim
      integer :: d, n, i, j, loop, k, ibest

      if (present(control)) then
         ctl = control
      else
         ctl = ps_control()
      end if

      d = size(xmin)
      n = ctl%n
      call validate_inputs(d, xmin, xmax, vmax, ctl)
      call rng%seed(ctl%seed)

      allocate(x(d,n), v(d,n), pbest(d,n), fx(n), fpbest(n), gbest(d))
      if (ctl%store_history) then
         allocate(result%best_history(ctl%max_loop), result%mean_history(ctl%max_loop))
      else
         allocate(result%best_history(0), result%mean_history(0))
      end if

      ! The original R code uses xmin[1], xmax[1] for every coordinate.
      do j = 1, d
         do i = 1, n
            if (ctl%legacy_initial_bounds) then
               x(j,i) = xmin(1) + (xmax(1) - xmin(1))*rng%uniform()
            else
               x(j,i) = xmin(j) + (xmax(j) - xmin(j))*rng%uniform()
            end if
         end do
      end do

      call evaluate_population(fun, x, fx, result%nfe)
      pbest = x
      fpbest = fx
      ibest = maxloc(fx, dim=1)
      gbest = x(:,ibest)
      gbest_value = fx(ibest)

      ! R's runif(n*d, min=-vmax, max=vmax) recycles vmax across the
      ! column-major linear draw stream. Preserve that behavior optionally.
      k = 0
      do j = 1, d
         do i = 1, n
            k = k + 1
            if (ctl%legacy_velocity_initialization) then
               vlim = vmax(1 + mod(k-1, d))
            else
               vlim = vmax(j)
            end if
            v(j,i) = -vlim + 2.0_dp*vlim*rng%uniform()
         end do
      end do

      do loop = 1, ctl%max_loop
         ! In the R implementation FUN(x) and FUN(pbest) are repeatedly
         ! reevaluated. With a pure pointwise objective, caching pbest values
         ! is mathematically equivalent and avoids needless evaluations.
         call evaluate_population(fun, x, fx, result%nfe)
         if (ctl%store_history) result%mean_history(loop) = sum(fx)/real(n,dp)

         do i = 1, n
            if (fx(i) > fpbest(i)) then
               pbest(:,i) = x(:,i)
               fpbest(i) = fx(i)
            end if
         end do

         ibest = maxloc(fpbest, dim=1)
         candidate_value = fpbest(ibest)
         if (candidate_value > gbest_value) then
            gbest = pbest(:,ibest)
            gbest_value = candidate_value
         end if
         if (ctl%store_history) result%best_history(loop) = gbest_value

         do i = 1, n
            do j = 1, d
               r1 = rng%uniform()
               r2 = rng%uniform()
               v(j,i) = ctl%w*v(j,i) + ctl%c1*r1*(pbest(j,i)-x(j,i)) &
                        + ctl%c2*r2*(gbest(j)-x(j,i))

               if (v(j,i) > vmax(j) .or. v(j,i) < -vmax(j)) then
                  if (ctl%legacy_velocity_clip) then
                     ! Literal psoptim 1.0 behavior: both signs become +vmax.
                     v(j,i) = vmax(j)
                  else
                     v(j,i) = max(-vmax(j), min(vmax(j), v(j,i)))
                  end if
               end if

               x_prev = x(j,i)
               x(j,i) = x(j,i) + v(j,i)
               if (x(j,i) > xmax(j)) x(j,i) = x_prev
               if (x(j,i) < xmin(j)) x(j,i) = x_prev
            end do
         end do
      end do

      allocate(result%sol(d))
      result%sol = gbest
      result%val = gbest_value
      result%iterations = ctl%max_loop
      allocate(result%final_population(d,n), result%final_velocity(d,n), result%personal_best(d,n))
      result%final_population = x
      result%final_velocity = v
      result%personal_best = pbest
   end subroutine ps_optimize

   subroutine evaluate_population(fun, x, fx, nfe)
      procedure(ps_scalar_objective) :: fun
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: fx(:)
      integer(int64), intent(inout) :: nfe
      integer :: i

      do i = 1, size(x,2)
         fx(i) = fun(x(:,i))
      end do
      nfe = nfe + int(size(x,2), int64)
   end subroutine evaluate_population

   subroutine validate_inputs(d, xmin, xmax, vmax, ctl)
      integer, intent(in) :: d
      real(dp), intent(in) :: xmin(:), xmax(:), vmax(:)
      type(ps_control), intent(in) :: ctl

      if (d < 1) error stop "ps_optimize: dimension must be positive"
      if (size(xmax) /= d) error stop "ps_optimize: xmin/xmax sizes differ"
      if (size(vmax) < d) error stop "ps_optimize: vmax must have at least d entries"
      if (any(xmax <= xmin)) error stop "ps_optimize: each xmax must exceed xmin"
      if (any(vmax(1:d) <= 0.0_dp)) error stop "ps_optimize: vmax must be positive"
      if (ctl%n < 1) error stop "ps_optimize: n must be positive"
      if (ctl%max_loop < 0) error stop "ps_optimize: max_loop must be nonnegative"
   end subroutine validate_inputs

end module psoptim
