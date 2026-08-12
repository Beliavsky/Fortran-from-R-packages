module flsss_types
  use flsss_kinds, only : dp, i8
  implicit none
  private

  type, public :: subset_solution
    integer, allocatable :: idx(:)
  end type subset_solution

  type, public :: subset_solutions
    type(subset_solution), allocatable :: sol(:)
    logical :: timed_out = .false.
    integer(i8) :: nodes = 0_i8
    integer(i8) :: pruned = 0_i8
    integer(i8) :: hash_lookups = 0_i8
    integer(i8) :: hash_candidates = 0_i8
    integer(i8) :: bound_states = 0_i8
    integer(i8) :: partitions_run = 0_i8
    integer(i8) :: tri_entries = 0_i8
    integer(i8) :: tri_lookups = 0_i8
    integer(i8) :: bound_updates = 0_i8
    integer(i8) :: mpat_splits = 0_i8
    integer :: packed_lanes = 0
    character(len=24) :: engine = 'dfs'
  contains
    procedure :: size => subset_solutions_size
  end type subset_solutions

  type, public :: real_bucket
    real(dp), allocatable :: value(:)
  end type real_bucket

  type, public :: multiset_solution
    type(subset_solution), allocatable :: bucket(:)
  end type multiset_solution

  type, public :: multiset_solutions
    type(multiset_solution), allocatable :: sol(:)
    logical :: timed_out = .false.
    integer(i8) :: nodes = 0_i8
  contains
    procedure :: size => multiset_solutions_size
  end type multiset_solutions

  type, public :: integerized_problem
    integer(i8), allocatable :: v(:,:)
    integer(i8), allocatable :: target(:)
    integer(i8), allocatable :: me(:)
    integer, allocatable :: ratio(:)
    integer :: compressed_dim = 0
  end type integerized_problem


  type, public :: integerized_search_result
    type(subset_solutions) :: solution
    type(integerized_problem) :: integerized
  end type integerized_search_result

  type, public :: knapsack_result
    integer, allocatable :: solution(:)
    real(dp), allocatable :: selection_costs(:)
    real(dp), allocatable :: budgets(:)
    real(dp) :: selection_profit = -huge(1.0_dp)
    real(dp) :: unconstrained_max_profit = 0.0_dp
    logical :: feasible = .false.
    logical :: timed_out = .false.
    integer(i8) :: nodes = 0_i8
    type(integerized_problem) :: integerized
  end type knapsack_result

  type, public :: knapsack_multi_result
    real(dp), allocatable :: max_value(:)
    type(subset_solution), allocatable :: selection(:)
    real(dp), allocatable :: lookup_table(:,:)
    logical :: timed_out = .false.
    integer(i8) :: nodes = 0_i8
  end type knapsack_multi_result

  type, public :: gap_result
    integer, allocatable :: assignment(:)
    real(dp), allocatable :: agent_cost(:)
    real(dp) :: total_profit_or_loss = 0.0_dp
    real(dp) :: unconstrained_max_profit = 0.0_dp
    logical :: feasible = .false.
    logical :: timed_out = .false.
    integer(i8) :: nodes = 0_i8
    integer(i8) :: bkp_solved = 0_i8
  end type gap_result

  type, public :: mflsss_object
    integer :: len = 0
    real(dp), allocatable :: v(:,:)
    real(dp), allocatable :: target(:), me(:)
    integer, allocatable :: lb(:), ub(:)
    integer :: dl = 0, du = 0
    integer :: prefix = 0
    integer :: prefix_lo = 0, prefix_hi = 0
  end type mflsss_object

  type, public :: mflsss_decomposition
    type(mflsss_object), allocatable :: object(:)
    type(subset_solutions) :: solutions_found
  end type mflsss_decomposition

  type, public :: arb_flsss_object
    integer :: len = 0
    character(len=:), allocatable :: v(:,:)
    character(len=:), allocatable :: target(:)
    integer :: prefix = 0
    integer :: prefix_lo = 0, prefix_hi = 0
  end type arb_flsss_object

  type, public :: arb_flsss_decomposition
    type(arb_flsss_object), allocatable :: object(:)
    type(subset_solutions) :: solutions_found
  end type arb_flsss_decomposition

  type, public :: ksum_table
    integer :: k = 0
    integer, allocatable :: index(:,:)
    character(len=:), allocatable :: sum(:,:)
    integer(i8), allocatable :: hash(:)
    integer, allocatable :: order(:)
    integer, allocatable :: frac_digits(:)
  end type ksum_table

contains

  integer function subset_solutions_size(this) result(n)
    class(subset_solutions), intent(in) :: this
    if (allocated(this%sol)) then
      n = size(this%sol)
    else
      n = 0
    end if
  end function subset_solutions_size

  integer function multiset_solutions_size(this) result(n)
    class(multiset_solutions), intent(in) :: this
    if (allocated(this%sol)) then
      n = size(this%sol)
    else
      n = 0
    end if
  end function multiset_solutions_size

end module flsss_types
