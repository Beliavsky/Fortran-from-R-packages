module trawl_optimize
  use trawl_kinds, only : dp
  use trawl_rng, only : runif_scalar
  use deoptim, only : de_control, de_result, deoptim_solve, de_success, i8
  implicit none
  private
  public :: differential_evolution

  abstract interface
    function objective_fn(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_fn
  end interface

contains

  subroutine differential_evolution(fun, lower, upper, best, bestval, itermax, npop, status, seed)
    procedure(objective_fn) :: fun
    real(dp), intent(in) :: lower(:), upper(:)
    real(dp), intent(out) :: best(:), bestval
    integer, intent(in), optional :: itermax, npop
    integer, intent(out), optional :: status
    integer(i8), intent(in), optional :: seed

    type(de_control) :: control
    type(de_result) :: result
    integer :: d, it

    d = size(lower)
    if (size(upper) /= d .or. size(best) /= d .or. d < 1) then
      best = 0.0_dp
      bestval = huge(1.0_dp)
      if (present(status)) status = 1
      return
    end if

    ! Match the DEoptim.control() defaults used by upstream trawl, except that
    ! trawl explicitly requests itermax=1000 and trace=FALSE.
    it = 1000
    if (present(itermax)) it = max(1, itermax)
    control%itermax = it
    control%strategy = 2
    control%cr = 0.5_dp
    control%f = 0.8_dp
    control%bs = .false.
    control%trace = 0
    control%np = 10 * d
    if (present(npop)) control%np = max(4, npop)
    control%p = 0.2_dp
    control%c = 0.0_dp
    control%reltol = sqrt(epsilon(1.0_dp))
    control%steptol = it
    control%storepopfrom = it + 1
    control%storepopfreq = 1
    if (present(seed)) then
      control%seed = seed
    else
      ! Upstream DEoptim draws from R's global RNG.  Consume one draw from the
      ! trawl RNG to seed the standalone DEoptim RNG so set_trawl_seed() also
      ! makes optimizer-backed fits reproducible.
      control%seed = 1_i8 + int(runif_scalar() * 2147483645.0_dp, i8)
    end if

    call deoptim_solve(fun, lower, upper, result, control=control)

    if (result%status == de_success .and. allocated(result%bestmem)) then
      best = result%bestmem
      bestval = result%bestval
    else
      best = 0.0_dp
      bestval = huge(1.0_dp)
    end if

    if (present(status)) then
      if (result%status == de_success .and. allocated(result%bestmem)) then
        status = 0
      else
        status = result%status
        if (status == 0) status = 1
      end if
    end if
  end subroutine differential_evolution

end module trawl_optimize
