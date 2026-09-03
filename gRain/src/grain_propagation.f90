module grain_propagation
  use r_kinds, only : dp
  use grbase_types, only : table_t
  use grbase_tables, only : valid_table, table_margin, table_multiply, table_divide_zero
  use grain_types, only : grain_network_t
  use grain_compile, only : grain_ok, grain_invalid_input, grain_impossible_evidence, grain_not_compiled
  implicit none
  private

  public :: propagate_ls
  public :: probability_of_evidence

contains

  pure subroutine propagate_ls(network, status)
    type(grain_network_t), intent(inout) :: network !! Compiled network whose temporary clique potentials are calibrated in place.
    integer, intent(out) :: status !! Zero on success; reports uncompiled or impossible-evidence states otherwise.
    type(table_t) :: separator_potential
    real(dp) :: normalization
    real(dp) :: z
    integer :: child
    integer :: i
    integer :: parent
    integer :: nc

    if (.not. network%compiled) then
      status = grain_not_compiled
      return
    end if
    if (.not. allocated(network%pot_temp)) then
      status = grain_invalid_input
      return
    end if
    nc = network%rip%cliques%count
    if (size(network%pot_temp) /= nc .or. nc < 1) then
      status = grain_invalid_input
      return
    end if

    network%pot_equi = network%pot_temp

    do i = nc, 2, -1
      parent = network%rip%parent(i)
      if (parent <= 0 .or. parent >= i) then
        status = grain_invalid_input
        return
      end if
      if (size(network%rip%separators%set(i)%value) > 0) then
        separator_potential = table_margin(network%pot_equi(i), &
                                           network%rip%separators%set(i)%value)
        if (.not. valid_table(separator_potential)) then
          status = grain_invalid_input
          return
        end if
        network%pot_equi(i) = table_divide_zero(network%pot_equi(i), separator_potential)
        network%pot_equi(parent) = table_multiply(network%pot_equi(parent), separator_potential)
      else
        z = sum(network%pot_equi(i)%value)
        if (z <= 0.0_dp) then
          network%p_evidence = 0.0_dp
          network%propagated = .false.
          status = grain_impossible_evidence
          return
        end if
        network%pot_equi(1)%value = network%pot_equi(1)%value * z
        network%pot_equi(i)%value = network%pot_equi(i)%value / z
      end if
    end do

    normalization = sum(network%pot_equi(1)%value)
    network%p_evidence = normalization
    if (normalization <= 0.0_dp) then
      network%propagated = .false.
      status = grain_impossible_evidence
      return
    end if
    network%pot_equi(1)%value = network%pot_equi(1)%value / normalization

    do i = 1, nc
      do child = i + 1, nc
        if (network%rip%parent(child) /= i) cycle
        if (size(network%rip%separators%set(child)%value) == 0) cycle
        separator_potential = table_margin(network%pot_equi(i), &
                                           network%rip%separators%set(child)%value)
        network%pot_equi(child) = table_multiply(network%pot_equi(child), separator_potential)
        if (.not. valid_table(network%pot_equi(child))) then
          status = grain_invalid_input
          return
        end if
      end do
    end do

    network%propagated = .true.
    status = grain_ok
  end subroutine propagate_ls

  function probability_of_evidence(network, status) result(probability)
    type(grain_network_t), intent(in) :: network !! Network whose current active evidence probability is requested.
    integer, intent(out) :: status !! Zero if the stored probability is valid, or `grain_not_compiled` otherwise.
    real(dp) :: probability

    if (.not. network%compiled) then
      probability = 0.0_dp
      status = grain_not_compiled
      return
    end if
    probability = network%p_evidence
    status = grain_ok
  end function probability_of_evidence

end module grain_propagation
