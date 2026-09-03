module grain_types
  use r_kinds, only : dp
  use grbase_types, only : table_t, rip_order_t
  implicit none
  private

  type, public :: grain_network_t
    integer :: n_nodes = 0
    integer, allocatable :: cardinality(:)
    integer, allocatable :: dag(:, :)
    integer, allocatable :: ug(:, :)
    type(rip_order_t) :: rip
    type(table_t), allocatable :: cpt(:)
    type(table_t), allocatable :: pot_orig(:)
    type(table_t), allocatable :: pot_temp(:)
    type(table_t), allocatable :: pot_equi(:)
    type(table_t), allocatable :: evidence(:)
    logical, allocatable :: evidence_active(:)
    logical, allocatable :: evidence_hard(:)
    integer, allocatable :: evidence_state(:)
    real(dp) :: p_evidence = 1.0_dp
    logical :: compiled = .false.
    logical :: propagated = .false.
  end type grain_network_t

end module grain_types
