module grbase_types
  use r_kinds, only : dp
  implicit none
  private

  type, public :: integer_set_t
    integer, allocatable :: value(:)
  end type integer_set_t

  type, public :: set_list_t
    type(integer_set_t), allocatable :: set(:)
    integer :: count = 0
  end type set_list_t

  type, public :: table_t
    integer, allocatable :: var(:)
    integer, allocatable :: dim(:)
    real(dp), allocatable :: value(:)
  end type table_t

  type, public :: rip_order_t
    type(set_list_t) :: cliques
    type(set_list_t) :: separators
    integer, allocatable :: parent(:)
    integer, allocatable :: child(:)
    integer, allocatable :: host(:)
    integer, allocatable :: order(:)
  end type rip_order_t

end module grbase_types
