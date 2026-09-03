module grain_update
  use r_kinds, only : dp
  use grain_types, only : grain_network_t
  use grain_cpt, only : replace_cpt_values
  use grain_compile, only : rebuild_potentials, grain_ok, grain_invalid_input, grain_not_compiled
  use grain_propagation, only : propagate_ls
  implicit none
  private

  public :: replace_network_cpt

contains

  pure subroutine replace_network_cpt(network, node, values, status, propagate)
    type(grain_network_t), intent(inout) :: network !! CPT network updated without changing its triangulation.
    integer, value :: node !! One-based child-node label indexing the CPT to replace.
    real(dp), intent(in) :: values(:) !! Replacement CPT weights in the existing table's flattened column-major order.
    integer, intent(out) :: status !! Zero on success; nonzero for an invalid node, length, or network state.
    logical, value, optional :: propagate !! If true, recalibrate after rebuilding clique potentials; defaults to false.
    logical :: do_propagate
    logical :: ok

    if (.not. network%compiled) then
      status = grain_not_compiled
      return
    end if
    if (.not. allocated(network%cpt)) then
      status = grain_invalid_input
      return
    end if
    if (node < 1 .or. node > size(network%cpt)) then
      status = grain_invalid_input
      return
    end if
    call replace_cpt_values(network%cpt(node), values, ok)
    if (.not. ok) then
      status = grain_invalid_input
      return
    end if
    call rebuild_potentials(network, status)
    if (status /= grain_ok) return
    do_propagate = .false.
    if (present(propagate)) do_propagate = propagate
    if (do_propagate) call propagate_ls(network, status)
  end subroutine replace_network_cpt

end module grain_update
