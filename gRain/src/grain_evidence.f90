module grain_evidence
  use r_kinds, only : dp
  use grbase_types, only : table_t
  use grbase_tables, only : make_table, valid_table
  use grain_types, only : grain_network_t
  use grain_compile, only : rebuild_potentials, grain_ok, grain_invalid_input, grain_not_compiled
  use grain_propagation, only : propagate_ls
  implicit none
  private

  public :: set_hard_evidence
  public :: set_soft_evidence
  public :: clear_evidence
  public :: absorb_evidence
  public :: has_evidence

contains

  pure subroutine set_hard_evidence(network, node, state, status, propagate, overwrite)
    type(grain_network_t), intent(inout) :: network !! Compiled network receiving a one-state indicator likelihood factor.
    integer, value :: node !! One-based node label on which hard evidence is asserted.
    integer, value :: state !! One-based state label that receives unit weight; all other states receive zero.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid node/state or an uncompiled network.
    logical, value, optional :: propagate !! If true, calibrate immediately after changing evidence; defaults to true.
    logical, value, optional :: overwrite !! Whether to replace prior node evidence; default false.
    real(dp), allocatable :: weights(:)
    logical :: do_propagate
    logical :: do_overwrite

    if (.not. network%compiled) then
      status = grain_not_compiled
      return
    end if
    if (node < 1 .or. node > network%n_nodes) then
      status = grain_invalid_input
      return
    end if
    if (state < 1 .or. state > network%cardinality(node)) then
      status = grain_invalid_input
      return
    end if
    do_overwrite = .false.
    if (present(overwrite)) do_overwrite = overwrite
    if (network%evidence_active(node) .and. .not. do_overwrite) then
      status = grain_ok
      return
    end if

    allocate(weights(network%cardinality(node)))
    weights = 0.0_dp
    weights(state) = 1.0_dp
    network%evidence(node) = make_table([node], [network%cardinality(node)], weights)
    network%evidence_active(node) = .true.
    network%evidence_hard(node) = .true.
    network%evidence_state(node) = state
    call rebuild_potentials(network, status)
    if (status /= grain_ok) return

    do_propagate = .true.
    if (present(propagate)) do_propagate = propagate
    if (do_propagate) call propagate_ls(network, status)
  end subroutine set_hard_evidence

  pure subroutine set_soft_evidence(network, node, weights, status, propagate, overwrite)
    type(grain_network_t), intent(inout) :: network !! Compiled network receiving a unary likelihood/virtual-evidence factor.
    integer, value :: node !! One-based node label on which the likelihood evidence is defined.
    real(dp), intent(in) :: weights(:) !! Nonnegative state likelihood weights; they need not sum to one but cannot all be zero.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid weights/node or an uncompiled network.
    logical, value, optional :: propagate !! If true, calibrate immediately after changing evidence; defaults to true.
    logical, value, optional :: overwrite !! If true, replace prior evidence on this node; defaults to false.
    logical :: do_propagate
    logical :: do_overwrite

    if (.not. network%compiled) then
      status = grain_not_compiled
      return
    end if
    if (node < 1 .or. node > network%n_nodes) then
      status = grain_invalid_input
      return
    end if
    if (size(weights) /= network%cardinality(node)) then
      status = grain_invalid_input
      return
    end if
    if (any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) then
      status = grain_invalid_input
      return
    end if
    do_overwrite = .false.
    if (present(overwrite)) do_overwrite = overwrite
    if (network%evidence_active(node) .and. .not. do_overwrite) then
      status = grain_ok
      return
    end if

    network%evidence(node) = make_table([node], [network%cardinality(node)], weights)
    if (.not. valid_table(network%evidence(node))) then
      status = grain_invalid_input
      return
    end if
    network%evidence_active(node) = .true.
    network%evidence_hard(node) = .false.
    network%evidence_state(node) = 0
    call rebuild_potentials(network, status)
    if (status /= grain_ok) return

    do_propagate = .true.
    if (present(propagate)) do_propagate = propagate
    if (do_propagate) call propagate_ls(network, status)
  end subroutine set_soft_evidence

  pure subroutine clear_evidence(network, status, node, propagate)
    type(grain_network_t), intent(inout) :: network !! Network from which all evidence or one node's evidence is removed.
    integer, intent(out) :: status !! Zero on success or `grain_not_compiled` for an uncompiled network.
    integer, value, optional :: node !! Optional one-based node whose evidence alone is cleared; absent means clear all evidence.
    logical, value, optional :: propagate !! If true, calibrate after clearing evidence; defaults to true.
    logical :: do_propagate

    if (.not. network%compiled) then
      status = grain_not_compiled
      return
    end if
    if (present(node)) then
      if (node < 1 .or. node > network%n_nodes) then
        status = grain_invalid_input
        return
      end if
      network%evidence_active(node) = .false.
      network%evidence_hard(node) = .false.
      network%evidence_state(node) = 0
      call clear_table(network%evidence(node))
    else
      network%evidence_active = .false.
      network%evidence_hard = .false.
      network%evidence_state = 0
      call clear_all_evidence_tables(network)
    end if

    call rebuild_potentials(network, status)
    if (status /= grain_ok) return
    do_propagate = .true.
    if (present(propagate)) do_propagate = propagate
    if (do_propagate) call propagate_ls(network, status)
  end subroutine clear_evidence

  pure subroutine absorb_evidence(network, status, propagate)
    type(grain_network_t), intent(inout) :: network !! Network whose evidence-weighted potentials become the baseline.
    integer, intent(out) :: status !! Zero on success or `grain_not_compiled` when the network has not been compiled.
    logical, value, optional :: propagate !! If true, calibrate the absorbed network immediately; defaults to true.
    logical :: do_propagate

    if (.not. network%compiled) then
      status = grain_not_compiled
      return
    end if
    network%pot_orig = network%pot_temp
    network%evidence_active = .false.
    network%evidence_hard = .false.
    network%evidence_state = 0
    call clear_all_evidence_tables(network)
    network%pot_temp = network%pot_orig
    network%pot_equi = network%pot_orig
    network%propagated = .false.
    network%p_evidence = 1.0_dp

    do_propagate = .true.
    if (present(propagate)) do_propagate = propagate
    if (do_propagate) then
      call propagate_ls(network, status)
    else
      status = grain_ok
    end if
  end subroutine absorb_evidence

  pure logical function has_evidence(network, node) result(active)
    type(grain_network_t), intent(in) :: network !! Network whose evidence bookkeeping is inspected.
    integer, value, optional :: node !! Optional node label; absent asks whether any evidence is active.

    active = .false.
    if (.not. allocated(network%evidence_active)) return
    if (present(node)) then
      if (node < 1 .or. node > size(network%evidence_active)) return
      active = network%evidence_active(node)
    else
      active = any(network%evidence_active)
    end if
  end function has_evidence

  pure subroutine clear_all_evidence_tables(network)
    type(grain_network_t), intent(inout) :: network !! Network whose allocated unary evidence tables are deallocated.
    integer :: i

    do i = 1, size(network%evidence)
      call clear_table(network%evidence(i))
    end do
  end subroutine clear_all_evidence_tables

  pure subroutine clear_table(tab)
    type(table_t), intent(inout) :: tab !! Table object whose allocatable metadata and values are released.

    if (allocated(tab%var)) deallocate(tab%var)
    if (allocated(tab%dim)) deallocate(tab%dim)
    if (allocated(tab%value)) deallocate(tab%value)
  end subroutine clear_table

end module grain_evidence
