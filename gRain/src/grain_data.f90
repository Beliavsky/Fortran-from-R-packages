module grain_data
  use r_kinds, only : dp
  use grbase_types, only : table_t, rip_order_t
  use grbase_arrays, only : cell_to_entry
  use grbase_tables, only : make_table, valid_table, table_normalize_all, table_normalize_first
  use grbase_tables, only : table_margin, table_multiply, table_divide_zero
  use grbase_graphs, only : is_dag, parents_of
  use grain_compile, only : grain_ok, grain_invalid_input, grain_inconsistent_cardinality, grain_not_dag
  implicit none
  private

  public :: empirical_table
  public :: estimate_cpt
  public :: estimate_cpts_from_data
  public :: estimate_clique_marginals
  public :: marginals_to_potentials
  public :: potentials_to_marginals

contains

  pure function empirical_table(data, vars, cardinality, smooth, normalize) result(tab)
    integer, intent(in) :: data(:, :) !! Case-by-node observations with one-based integer states.
    integer, intent(in) :: vars(:) !! Distinct node labels whose joint empirical table is requested, in output dimension order.
    integer, intent(in) :: cardinality(:) !! State count for every data column/node.
    real(dp), value, optional :: smooth !! Nonnegative pseudocount added to every cell before optional normalization; default zero.
    logical, value, optional :: normalize !! If true, normalize the whole table to unit mass; defaults to false.
    type(table_t) :: tab
    integer, allocatable :: cell(:)
    integer, allocatable :: dim(:)
    real(dp), allocatable :: counts(:)
    real(dp) :: pseudocount
    integer :: entry
    integer :: i
    integer :: j
    logical :: do_normalize
    logical :: row_valid

    if (size(vars) < 1) return
    if (size(cardinality) /= size(data, 2)) return
    if (any(cardinality <= 0)) return
    if (any(vars < 1) .or. any(vars > size(data, 2))) return
    if (has_duplicates(vars)) return
    pseudocount = 0.0_dp
    if (present(smooth)) pseudocount = smooth
    if (pseudocount < 0.0_dp) return

    dim = cardinality(vars)
    allocate(counts(product(dim)), cell(size(vars)))
    counts = pseudocount
    do i = 1, size(data, 1)
      row_valid = .true.
      do j = 1, size(vars)
        cell(j) = data(i, vars(j))
        if (cell(j) < 1 .or. cell(j) > dim(j)) row_valid = .false.
      end do
      if (.not. row_valid) cycle
      entry = cell_to_entry(cell, dim)
      counts(entry) = counts(entry) + 1.0_dp
    end do
    tab = make_table(vars, dim, counts)
    do_normalize = .false.
    if (present(normalize)) do_normalize = normalize
    if (do_normalize .and. valid_table(tab)) tab = table_normalize_all(tab)
  end function empirical_table

  pure function estimate_cpt(data, child, parents, cardinality, smooth) result(cpt)
    integer, intent(in) :: data(:, :) !! Observation matrix with one-based categorical states and one column per node.
    integer, value :: child !! One-based child-node label for the conditional distribution.
    integer, intent(in) :: parents(:) !! Parent-node labels in the desired conditioning order.
    integer, intent(in) :: cardinality(:) !! State counts for all nodes/data columns.
    real(dp), value, optional :: smooth !! Nonnegative cell pseudocount added before conditional normalization; default zero.
    type(table_t) :: cpt
    type(table_t) :: counts
    integer, allocatable :: vars(:)
    real(dp) :: pseudocount

    if (child < 1 .or. child > size(data, 2)) return
    allocate(vars(1 + size(parents)))
    vars(1) = child
    if (size(parents) > 0) vars(2:) = parents
    pseudocount = 0.0_dp
    if (present(smooth)) pseudocount = smooth
    counts = empirical_table(data, vars, cardinality, pseudocount, normalize=.false.)
    if (.not. valid_table(counts)) return
    cpt = table_normalize_first(counts)
  end function estimate_cpt

  pure subroutine estimate_cpts_from_data(data, dag, cardinality, cpts, status, smooth)
    integer, intent(in) :: data(:, :) !! Observation matrix with rows as cases and columns corresponding to node labels 1..n.
    integer, intent(in) :: dag(:, :) !! Parent-to-child zero-one adjacency matrix for a directed acyclic graph.
    integer, intent(in) :: cardinality(:) !! State counts for all nodes; length must equal the DAG order and data columns.
    type(table_t), allocatable, intent(out) :: cpts(:) !! Estimated CPTs indexed by child node.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, cardinalities, or a cyclic graph.
    real(dp), value, optional :: smooth !! Nonnegative cell pseudocount applied to every estimated CPT; default zero.
    integer, allocatable :: parents(:)
    real(dp) :: pseudocount
    integer :: node
    integer :: n

    status = grain_invalid_input
    n = size(dag, 1)
    if (size(dag, 2) /= n .or. size(data, 2) /= n .or. size(cardinality) /= n) return
    if (any(cardinality <= 0)) then
      status = grain_inconsistent_cardinality
      return
    end if
    if (.not. is_dag(dag)) then
      status = grain_not_dag
      return
    end if
    pseudocount = 0.0_dp
    if (present(smooth)) pseudocount = smooth
    if (pseudocount < 0.0_dp) return

    allocate(cpts(n))
    do node = 1, n
      parents = parents_of(dag, [node])
      cpts(node) = estimate_cpt(data, node, parents, cardinality, pseudocount)
      if (.not. valid_table(cpts(node))) return
    end do
    status = grain_ok
  end subroutine estimate_cpts_from_data

  pure subroutine estimate_clique_marginals(data, rip_order, cardinality, marginals, status, smooth)
    integer, intent(in) :: data(:, :) !! Observation matrix with rows as cases and node columns using one-based states.
    type(rip_order_t), intent(in) :: rip_order !! Running-intersection order whose clique empirical marginals are estimated.
    integer, intent(in) :: cardinality(:) !! State counts for all nodes/data columns.
    type(table_t), allocatable, intent(out) :: marginals(:) !! Unit-mass empirical clique tables in RIP order.
    integer, intent(out) :: status !! Zero when every clique table is estimated successfully.
    real(dp), value, optional :: smooth !! Nonnegative pseudocount added to each clique cell; default zero.
    real(dp) :: pseudocount
    integer :: i

    status = grain_invalid_input
    if (size(data, 2) /= size(cardinality)) return
    pseudocount = 0.0_dp
    if (present(smooth)) pseudocount = smooth
    if (pseudocount < 0.0_dp) return
    allocate(marginals(rip_order%cliques%count))
    do i = 1, rip_order%cliques%count
      marginals(i) = empirical_table(data, rip_order%cliques%set(i)%value, cardinality, &
                                     pseudocount, normalize=.true.)
      if (.not. valid_table(marginals(i))) return
    end do
    status = grain_ok
  end subroutine estimate_clique_marginals

  pure subroutine marginals_to_potentials(marginals, rip_order, potentials, status)
    type(table_t), intent(in) :: marginals(:) !! Clique marginal tables in RIP order.
    type(rip_order_t), intent(in) :: rip_order !! Running-intersection order providing separator domains.
    type(table_t), allocatable, intent(out) :: potentials(:) !! Root marginal then clique-given-separator potentials.
    integer, intent(out) :: status !! Zero when all clique and separator operations are valid.
    type(table_t) :: separator
    integer :: i

    status = grain_invalid_input
    if (size(marginals) /= rip_order%cliques%count) return
    allocate(potentials(size(marginals)))
    do i = 1, size(marginals)
      if (.not. valid_table(marginals(i))) return
      if (i == 1 .or. size(rip_order%separators%set(i)%value) == 0) then
        potentials(i) = marginals(i)
      else
        separator = table_margin(marginals(i), rip_order%separators%set(i)%value)
        potentials(i) = table_divide_zero(marginals(i), separator)
      end if
      if (.not. valid_table(potentials(i))) return
    end do
    status = grain_ok
  end subroutine marginals_to_potentials

  pure subroutine potentials_to_marginals(potentials, rip_order, marginals, status)
    type(table_t), intent(in) :: potentials(:) !! Root marginal and clique-given-separator potentials in RIP order.
    type(rip_order_t), intent(in) :: rip_order !! Running-intersection order linking each nonroot clique to its parent.
    type(table_t), allocatable, intent(out) :: marginals(:) !! Reconstructed clique marginals in RIP order.
    integer, intent(out) :: status !! Zero when all parent/separator metadata and table operations are valid.
    type(table_t) :: separator
    integer :: i
    integer :: parent

    status = grain_invalid_input
    if (size(potentials) /= rip_order%cliques%count .or. size(potentials) < 1) return
    allocate(marginals(size(potentials)))
    marginals(1) = potentials(1)
    if (.not. valid_table(marginals(1))) return
    do i = 2, size(potentials)
      if (.not. valid_table(potentials(i))) return
      parent = rip_order%parent(i)
      if (parent < 1 .or. parent >= i) return
      if (size(rip_order%separators%set(i)%value) == 0) then
        marginals(i) = potentials(i)
      else
        separator = table_margin(marginals(parent), rip_order%separators%set(i)%value)
        marginals(i) = table_multiply(potentials(i), separator)
      end if
      if (.not. valid_table(marginals(i))) return
    end do
    status = grain_ok
  end subroutine potentials_to_marginals

  pure logical function has_duplicates(values) result(duplicate)
    integer, intent(in) :: values(:) !! Integer labels checked for repeated entries.
    integer :: i
    integer :: j

    duplicate = .false.
    do i = 1, size(values) - 1
      do j = i + 1, size(values)
        if (values(i) == values(j)) then
          duplicate = .true.
          return
        end if
      end do
    end do
  end function has_duplicates

end module grain_data
