module grain_query
  use r_kinds, only : dp
  use grbase_types, only : table_t
  use grbase_tables, only : valid_table, table_margin, table_multiply, table_divide_zero
  use grbase_sets, only : is_subset_of
  use grain_types, only : grain_network_t
  use grain_compile, only : grain_ok, grain_invalid_input, grain_not_compiled
  use grain_propagation, only : propagate_ls
  implicit none
  private

  public :: query_marginal
  public :: query_all_marginals
  public :: query_joint
  public :: query_conditional
  public :: reconstructed_joint

contains

  pure subroutine query_marginal(network, node, marginal, status)
    type(grain_network_t), intent(inout) :: network !! Network queried and propagated automatically when needed.
    integer, value :: node !! One-based node label whose posterior marginal distribution is requested.
    type(table_t), intent(out) :: marginal !! One-dimensional normalized posterior table for `node`.
    integer, intent(out) :: status !! Zero on success; nonzero for an invalid node or network state.
    integer :: host
    real(dp) :: total

    call ensure_propagated(network, status)
    if (status /= grain_ok) return
    if (node < 1 .or. node > network%n_nodes) then
      status = grain_invalid_input
      return
    end if
    host = network%rip%host(node)
    if (host < 1 .or. host > network%rip%cliques%count) then
      host = first_host(network, [node])
    end if
    if (host == 0) then
      status = grain_invalid_input
      return
    end if
    marginal = table_margin(network%pot_equi(host), [node])
    if (.not. valid_table(marginal)) then
      status = grain_invalid_input
      return
    end if
    total = sum(marginal%value)
    if (total > 0.0_dp) marginal%value = marginal%value / total
    status = grain_ok
  end subroutine query_marginal

  pure subroutine query_all_marginals(network, marginals, status)
    type(grain_network_t), intent(inout) :: network !! Network queried for every one-node posterior marginal.
    type(table_t), allocatable, intent(out) :: marginals(:) !! Array of one-dimensional tables indexed by node label.
    integer, intent(out) :: status !! Zero if all node marginals were computed successfully.
    integer :: i

    allocate(marginals(network%n_nodes))
    do i = 1, network%n_nodes
      call query_marginal(network, i, marginals(i), status)
      if (status /= grain_ok) return
    end do
  end subroutine query_all_marginals

  subroutine query_joint(network, nodes, joint, status, normalize)
    type(grain_network_t), intent(inout) :: network !! Network queried and propagated automatically when needed.
    integer, intent(in) :: nodes(:) !! Distinct one-based node labels whose joint posterior is requested in this output order.
    type(table_t), intent(out) :: joint !! Joint posterior table over `nodes`.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid or duplicate node labels.
    logical, value, optional :: normalize !! If false, return joint probability times evidence probability; defaults to true.
    type(table_t) :: full
    integer :: host
    logical :: do_normalize
    real(dp) :: total

    call ensure_propagated(network, status)
    if (status /= grain_ok) return
    if (.not. valid_node_set(network, nodes)) then
      status = grain_invalid_input
      return
    end if

    host = first_host(network, nodes)
    if (host > 0) then
      joint = table_margin(network%pot_equi(host), nodes)
    else
      full = reconstructed_joint(network, status)
      if (status /= grain_ok) return
      joint = table_margin(full, nodes)
    end if
    if (.not. valid_table(joint)) then
      status = grain_invalid_input
      return
    end if

    do_normalize = .true.
    if (present(normalize)) do_normalize = normalize
    if (do_normalize) then
      total = sum(joint%value)
      if (total > 0.0_dp) joint%value = joint%value / total
    else
      joint%value = joint%value * network%p_evidence
    end if
    status = grain_ok
  end subroutine query_joint

  subroutine query_conditional(network, nodes, conditional, status)
    type(grain_network_t), intent(inout) :: network !! Network queried for P(nodes(1) | nodes(2:)).
    integer, intent(in) :: nodes(:) !! Child node first followed by conditioning nodes; all labels must be distinct.
    type(table_t), intent(out) :: conditional !! Table normalized over its first dimension for each conditioning state.
    integer, intent(out) :: status !! Zero on success; invalid-input status when `nodes` is empty or malformed.
    type(table_t) :: joint
    type(table_t) :: conditioning

    if (size(nodes) < 1) then
      status = grain_invalid_input
      return
    end if
    call query_joint(network, nodes, joint, status, normalize=.true.)
    if (status /= grain_ok) return
    if (size(nodes) == 1) then
      conditional = joint
      return
    end if
    conditioning = table_margin(joint, nodes(2:))
    conditional = table_divide_zero(joint, conditioning)
    if (.not. valid_table(conditional)) then
      status = grain_invalid_input
      return
    end if
    status = grain_ok
  end subroutine query_conditional

  function reconstructed_joint(network, status) result(joint)
    type(grain_network_t), intent(in) :: network !! Propagated junction tree whose posterior joint is reconstructed.
    integer, intent(out) :: status !! Zero when the calibrated clique/separator factorization is valid.
    type(table_t) :: joint
    type(table_t) :: separator
    integer :: i

    if (.not. network%compiled .or. .not. network%propagated) then
      status = grain_not_compiled
      return
    end if
    if (network%rip%cliques%count < 1) then
      status = grain_invalid_input
      return
    end if
    joint = network%pot_equi(1)
    do i = 2, network%rip%cliques%count
      joint = table_multiply(joint, network%pot_equi(i))
      if (size(network%rip%separators%set(i)%value) > 0) then
        separator = table_margin(network%pot_equi(i), network%rip%separators%set(i)%value)
        joint = table_divide_zero(joint, separator)
      end if
      if (.not. valid_table(joint)) then
        status = grain_invalid_input
        return
      end if
    end do
    status = grain_ok
  end function reconstructed_joint

  pure subroutine ensure_propagated(network, status)
    type(grain_network_t), intent(inout) :: network !! Network calibrated if it is compiled but not yet propagated.
    integer, intent(out) :: status !! Zero after successful calibration, or an explanatory network status otherwise.

    if (.not. network%compiled) then
      status = grain_not_compiled
    else if (.not. network%propagated) then
      call propagate_ls(network, status)
    else
      status = grain_ok
    end if
  end subroutine ensure_propagated

  pure integer function first_host(network, nodes) result(host)
    type(grain_network_t), intent(in) :: network !! Compiled network whose ordered clique list is searched.
    integer, intent(in) :: nodes(:) !! Node set that must be contained in the returned clique.
    integer :: i

    host = 0
    do i = 1, network%rip%cliques%count
      if (is_subset_of(nodes, network%rip%cliques%set(i)%value)) then
        host = i
        return
      end if
    end do
  end function first_host

  pure logical function valid_node_set(network, nodes) result(ok)
    type(grain_network_t), intent(in) :: network !! Network defining the valid node-label range.
    integer, intent(in) :: nodes(:) !! Candidate list of distinct node labels.
    integer :: i
    integer :: j

    ok = size(nodes) > 0
    if (.not. ok) return
    ok = all(nodes >= 1 .and. nodes <= network%n_nodes)
    if (.not. ok) return
    do i = 1, size(nodes) - 1
      do j = i + 1, size(nodes)
        if (nodes(i) == nodes(j)) then
          ok = .false.
          return
        end if
      end do
    end do
  end function valid_node_set

end module grain_query
