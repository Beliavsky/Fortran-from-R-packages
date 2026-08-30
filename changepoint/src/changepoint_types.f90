! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_types
use r_kinds, only : dp
implicit none
private

integer, parameter, public :: cp_ok = 0
integer, parameter, public :: cp_invalid_argument = 1
integer, parameter, public :: cp_invalid_data = 2
integer, parameter, public :: cp_linalg_failure = 3

type, public :: changepoint_result
    integer :: status = cp_ok
    integer :: ncpts = 0
    real(dp) :: objective = 0.0_dp
    real(dp) :: unpenalized_cost = 0.0_dp
    integer, allocatable :: cpts(:)
end type changepoint_result

type, public :: amoc_result
    integer :: status = cp_ok
    integer :: cpt = 0
    logical :: changed = .false.
    real(dp) :: null_cost = 0.0_dp
    real(dp) :: alt_cost = 0.0_dp
    real(dp) :: test_statistic = 0.0_dp
end type amoc_result

type, public :: binseg_result
    integer :: status = cp_ok
    integer :: ncpts = 0
    integer, allocatable :: candidates(:)
    real(dp), allocatable :: scores(:)
    integer, allocatable :: cpts(:)
end type binseg_result

type, public :: segneigh_result
    integer :: status = cp_ok
    integer :: ncpts = 0
    real(dp), allocatable :: cost_by_ncpts(:)
    integer, allocatable :: cpts(:)
end type segneigh_result

type, public :: crops_solution
    real(dp) :: beta_start = 0.0_dp
    real(dp) :: beta_end = 0.0_dp
    real(dp) :: unpenalized_cost = 0.0_dp
    integer :: ncpts = 0
    integer, allocatable :: cpts(:)
end type crops_solution

end module changepoint_types
