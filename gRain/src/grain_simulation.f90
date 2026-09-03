module grain_simulation
  use, intrinsic :: iso_fortran_env, only : int64
  use r_kinds, only : dp
  use grbase_types, only : table_t
  use grbase_tables, only : valid_table, table_margin, table_divide_zero, table_slice
  use grbase_tables, only : table_sample_from_uniforms
  use grain_types, only : grain_network_t
  use grain_compile, only : grain_ok, grain_invalid_input
  use grain_propagation, only : propagate_ls
  implicit none
  private

  public :: simulate_from_uniforms
  public :: simulate_network

contains

  pure subroutine simulate_from_uniforms(network, uniforms, samples, status)
    type(grain_network_t), intent(inout) :: network !! Network sampled from its propagated posterior distribution.
    real(dp), intent(in) :: uniforms(:, :) !! Uniform draws with shape `(nsim,ncliques)`; values are clipped to [0,1].
    integer, allocatable, intent(out) :: samples(:, :) !! Simulated one-based state indices with shape `(nsim,n_nodes)`.
    integer, intent(out) :: status !! Zero on success; invalid-input status for shape/domain inconsistencies.
    type(table_t) :: conditional
    type(table_t) :: separator_marginal
    type(table_t) :: sliced
    integer, allocatable :: cells(:, :)
    integer, allocatable :: fixed_levels(:)
    integer :: c
    integer :: i
    integer :: j
    integer :: nsim
    integer :: v

    if (.not. network%propagated) then
      call propagate_ls(network, status)
      if (status /= grain_ok) return
    end if
    nsim = size(uniforms, 1)
    if (size(uniforms, 2) /= network%rip%cliques%count .or. nsim < 0) then
      status = grain_invalid_input
      return
    end if
    allocate(samples(nsim, network%n_nodes))
    samples = 0
    if (nsim == 0) then
      status = grain_ok
      return
    end if

    cells = table_sample_from_uniforms(network%pot_equi(1), uniforms(:, 1))
    if (size(cells, 1) /= nsim .or. size(cells, 2) /= size(network%pot_equi(1)%var)) then
      status = grain_invalid_input
      return
    end if
    do j = 1, size(network%pot_equi(1)%var)
      v = network%pot_equi(1)%var(j)
      samples(:, v) = cells(:, j)
    end do

    do c = 2, network%rip%cliques%count
      if (size(network%rip%separators%set(c)%value) > 0) then
        separator_marginal = table_margin(network%pot_equi(c), network%rip%separators%set(c)%value)
        conditional = table_divide_zero(network%pot_equi(c), separator_marginal)
      else
        conditional = network%pot_equi(c)
      end if
      if (.not. valid_table(conditional)) then
        status = grain_invalid_input
        return
      end if

      do i = 1, nsim
        if (size(network%rip%separators%set(c)%value) > 0) then
          allocate(fixed_levels(size(network%rip%separators%set(c)%value)))
          do j = 1, size(fixed_levels)
            fixed_levels(j) = samples(i, network%rip%separators%set(c)%value(j))
          end do
          sliced = table_slice(conditional, network%rip%separators%set(c)%value, fixed_levels)
          deallocate(fixed_levels)
        else
          sliced = conditional
        end if
        if (.not. valid_table(sliced)) then
          status = grain_invalid_input
          return
        end if
        cells = table_sample_from_uniforms(sliced, [uniforms(i, c)])
        if (size(cells, 1) /= 1 .or. size(cells, 2) /= size(sliced%var)) then
          status = grain_invalid_input
          return
        end if
        do j = 1, size(sliced%var)
          samples(i, sliced%var(j)) = cells(1, j)
        end do
      end do
    end do

    if (any(samples < 1)) then
      status = grain_invalid_input
      return
    end if
    status = grain_ok
  end subroutine simulate_from_uniforms

  pure subroutine simulate_network(network, nsim, samples, status, seed)
    type(grain_network_t), intent(inout) :: network !! Network sampled with the portable internal random stream.
    integer, value :: nsim !! Number of posterior cases to simulate; must be nonnegative.
    integer, allocatable, intent(out) :: samples(:, :) !! Simulated one-based state indices with one row per requested case.
    integer, intent(out) :: status !! Zero on success; invalid input for a negative sample count.
    integer, value, optional :: seed !! Positive integer seed for the portable Park-Miller stream; default 12345.
    real(dp), allocatable :: uniforms(:, :)
    integer(int64) :: state
    integer :: c
    integer :: i
    integer :: seed_value

    if (nsim < 0) then
      status = grain_invalid_input
      allocate(samples(0, max(0, network%n_nodes)))
      return
    end if
    seed_value = 12345
    if (present(seed)) seed_value = seed
    if (seed_value <= 0) seed_value = 12345
    state = int(mod(seed_value, 2147483646), int64) + 1_int64
    allocate(uniforms(nsim, network%rip%cliques%count))
    do c = 1, size(uniforms, 2)
      do i = 1, nsim
        state = modulo(16807_int64 * state, 2147483647_int64)
        uniforms(i, c) = real(state, dp) / 2147483647.0_dp
      end do
    end do
    call simulate_from_uniforms(network, uniforms, samples, status)
  end subroutine simulate_network

end module grain_simulation
