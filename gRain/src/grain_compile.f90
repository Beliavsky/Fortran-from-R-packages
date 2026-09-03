module grain_compile
  use r_kinds, only : dp
  use grbase_types, only : table_t, rip_order_t
  use grbase_tables, only : make_table, valid_table, table_multiply
  use grbase_graphs, only : is_dag, moralize_graph, minimal_triangulation
  use grbase_decompositions, only : rip_from_adjacency
  use grbase_sets, only : is_subset_of
  use grain_types, only : grain_network_t
  implicit none
  private

  integer, parameter, public :: grain_ok = 0
  integer, parameter, public :: grain_invalid_input = 1
  integer, parameter, public :: grain_inconsistent_cardinality = 2
  integer, parameter, public :: grain_not_dag = 3
  integer, parameter, public :: grain_impossible_evidence = 4
  integer, parameter, public :: grain_not_compiled = 5

  public :: compile_network
  public :: rebuild_potentials
  public :: initialize_network_from_potentials

contains

  subroutine compile_network(cpts, network, status, root)
    type(table_t), intent(in) :: cpts(:) !! One CPT per node; each CPT lists its child first and its parents after it.
    type(grain_network_t), intent(out) :: network !! Compiled Bayesian network with triangulation, RIP order, and clique potentials.
    integer, intent(out) :: status !! Zero on success; otherwise one of the public `grain_*` status constants.
    integer, intent(in), optional :: root(:) !! Optional node set forced complete before triangulation so one clique can contain it.
    integer, allocatable :: cardinality(:)
    integer, allocatable :: dag(:, :)
    integer, allocatable :: moral(:, :)
    integer, allocatable :: tug(:, :)
    type(table_t), allocatable :: ordered_cpt(:)
    integer :: child
    integer :: i
    integer :: j
    integer :: n

    status = grain_invalid_input
    n = size(cpts)
    if (n < 1) return
    allocate(cardinality(n), dag(n, n), ordered_cpt(n))
    cardinality = 0
    dag = 0

    do i = 1, n
      if (.not. valid_table(cpts(i))) return
      if (size(cpts(i)%var) < 1) return
      child = cpts(i)%var(1)
      if (child < 1 .or. child > n) return
      if (allocated(ordered_cpt(child)%var)) return
      ordered_cpt(child) = cpts(i)
      do j = 1, size(cpts(i)%var)
        if (cpts(i)%var(j) < 1 .or. cpts(i)%var(j) > n) return
        if (cardinality(cpts(i)%var(j)) == 0) then
          cardinality(cpts(i)%var(j)) = cpts(i)%dim(j)
        else if (cardinality(cpts(i)%var(j)) /= cpts(i)%dim(j)) then
          status = grain_inconsistent_cardinality
          return
        end if
      end do
      do j = 2, size(cpts(i)%var)
        if (cpts(i)%var(j) == child) return
        dag(cpts(i)%var(j), child) = 1
      end do
    end do

    if (any(cardinality == 0)) return
    do i = 1, n
      if (.not. allocated(ordered_cpt(i)%var)) return
      if (ordered_cpt(i)%dim(1) /= cardinality(i)) then
        status = grain_inconsistent_cardinality
        return
      end if
    end do
    if (.not. is_dag(dag)) then
      status = grain_not_dag
      return
    end if

    moral = moralize_graph(dag)
    if (present(root)) then
      if (any(root < 1) .or. any(root > n)) return
      call complete_vertex_set(moral, root)
    end if
    tug = minimal_triangulation(moral, real(cardinality, dp))
    if (size(tug, 1) /= n) return

    network%n_nodes = n
    network%cardinality = cardinality
    network%dag = dag
    network%ug = tug
    network%rip = rip_from_adjacency(tug)
    if (network%rip%cliques%count < 1) return
    network%cpt = ordered_cpt
    call allocate_evidence_state(network)
    call rebuild_potentials(network, status)
    if (status /= grain_ok) return
    network%compiled = .true.
    network%propagated = .false.
    network%p_evidence = 1.0_dp
  end subroutine compile_network

  pure subroutine rebuild_potentials(network, status)
    type(grain_network_t), intent(inout) :: network !! Network whose clique potentials are rebuilt from its current normalized CPTs.
    integer, intent(out) :: status !! Zero on success; nonzero if a CPT cannot be hosted by any clique or metadata are invalid.
    type(table_t) :: unity
    real(dp), allocatable :: ones(:)
    integer :: host
    integer :: i
    integer :: j
    integer :: nc

    status = grain_invalid_input
    if (.not. allocated(network%cpt)) return
    if (.not. allocated(network%cardinality)) return
    nc = network%rip%cliques%count
    if (nc < 1) return

    if (allocated(network%pot_orig)) deallocate(network%pot_orig)
    if (allocated(network%pot_temp)) deallocate(network%pot_temp)
    if (allocated(network%pot_equi)) deallocate(network%pot_equi)
    allocate(network%pot_orig(nc), network%pot_temp(nc), network%pot_equi(nc))

    do i = 1, nc
      allocate(ones(product(network%cardinality(network%rip%cliques%set(i)%value))))
      ones = 1.0_dp
      unity = make_table(network%rip%cliques%set(i)%value, &
                         network%cardinality(network%rip%cliques%set(i)%value), ones)
      deallocate(ones)
      if (.not. valid_table(unity)) return
      network%pot_orig(i) = unity
    end do

    do i = 1, size(network%cpt)
      host = 0
      do j = 1, nc
        if (is_subset_of(network%cpt(i)%var, network%rip%cliques%set(j)%value)) then
          host = j
          exit
        end if
      end do
      if (host == 0) return
      network%pot_orig(host) = table_multiply(network%pot_orig(host), network%cpt(i))
      if (.not. valid_table(network%pot_orig(host))) return
    end do

    network%pot_temp = network%pot_orig
    network%pot_equi = network%pot_orig
    call apply_active_evidence(network, status)
    if (status /= grain_ok) return
    network%propagated = .false.
    network%p_evidence = 1.0_dp
    status = grain_ok
  end subroutine rebuild_potentials

  pure subroutine initialize_network_from_potentials(potentials, rip_order, cardinality, network, status)
    type(table_t), intent(in) :: potentials(:) !! Clique-conditionals or clique potentials in the supplied RIP order.
    type(rip_order_t), intent(in) :: rip_order !! Valid running-intersection ordering describing the potential domains.
    integer, intent(in) :: cardinality(:) !! Per-node state counts for integer node labels 1..n.
    type(grain_network_t), intent(out) :: network !! Network initialized from clique potentials without DAG/CPT metadata.
    integer, intent(out) :: status !! Zero when all potential domains and cardinalities are consistent.
    integer :: i
    integer :: j

    status = grain_invalid_input
    if (size(potentials) /= rip_order%cliques%count) return
    if (size(cardinality) < 1 .or. any(cardinality <= 0)) return
    do i = 1, size(potentials)
      if (.not. valid_table(potentials(i))) return
      do j = 1, size(potentials(i)%var)
        if (potentials(i)%var(j) < 1 .or. potentials(i)%var(j) > size(cardinality)) return
        if (potentials(i)%dim(j) /= cardinality(potentials(i)%var(j))) then
          status = grain_inconsistent_cardinality
          return
        end if
      end do
    end do

    network%n_nodes = size(cardinality)
    network%cardinality = cardinality
    network%rip = rip_order
    allocate(network%dag(network%n_nodes, network%n_nodes))
    allocate(network%ug(network%n_nodes, network%n_nodes))
    network%dag = 0
    network%ug = 0
    network%pot_orig = potentials
    network%pot_temp = potentials
    network%pot_equi = potentials
    call allocate_evidence_state(network)
    network%compiled = .true.
    network%propagated = .false.
    network%p_evidence = 1.0_dp
    status = grain_ok
  end subroutine initialize_network_from_potentials

  pure subroutine complete_vertex_set(adj, vertices)
    integer, intent(inout) :: adj(:, :) !! Symmetric adjacency matrix modified by adding every edge among `vertices`.
    integer, intent(in) :: vertices(:) !! One-based node labels that must form a complete set.
    integer :: i
    integer :: j

    do i = 1, size(vertices) - 1
      do j = i + 1, size(vertices)
        if (vertices(i) == vertices(j)) cycle
        adj(vertices(i), vertices(j)) = 1
        adj(vertices(j), vertices(i)) = 1
      end do
    end do
  end subroutine complete_vertex_set

  pure subroutine allocate_evidence_state(network)
    type(grain_network_t), intent(inout) :: network !! Network receiving empty evidence bookkeeping arrays sized to its node count.
    integer :: i

    if (allocated(network%evidence)) deallocate(network%evidence)
    if (allocated(network%evidence_active)) deallocate(network%evidence_active)
    if (allocated(network%evidence_hard)) deallocate(network%evidence_hard)
    if (allocated(network%evidence_state)) deallocate(network%evidence_state)
    allocate(network%evidence(network%n_nodes))
    allocate(network%evidence_active(network%n_nodes))
    allocate(network%evidence_hard(network%n_nodes))
    allocate(network%evidence_state(network%n_nodes))
    network%evidence_active = .false.
    network%evidence_hard = .false.
    network%evidence_state = 0
    do i = 1, network%n_nodes
      if (allocated(network%evidence(i)%var)) deallocate(network%evidence(i)%var)
      if (allocated(network%evidence(i)%dim)) deallocate(network%evidence(i)%dim)
      if (allocated(network%evidence(i)%value)) deallocate(network%evidence(i)%value)
    end do
  end subroutine allocate_evidence_state

  pure subroutine apply_active_evidence(network, status)
    type(grain_network_t), intent(inout) :: network !! Network receiving active evidence factors in its temporary potentials.
    integer, intent(out) :: status !! Zero on success; invalid input if an active evidence factor has no host clique.
    integer :: host
    integer :: i
    integer :: j

    status = grain_ok
    if (.not. allocated(network%evidence_active)) return
    do i = 1, network%n_nodes
      if (.not. network%evidence_active(i)) cycle
      host = 0
      do j = 1, network%rip%cliques%count
        if (any(network%rip%cliques%set(j)%value == i)) then
          host = j
          exit
        end if
      end do
      if (host == 0) then
        status = grain_invalid_input
        return
      end if
      network%pot_temp(host) = table_multiply(network%pot_temp(host), network%evidence(i))
      if (.not. valid_table(network%pot_temp(host))) then
        status = grain_invalid_input
        return
      end if
    end do
  end subroutine apply_active_evidence

end module grain_compile
