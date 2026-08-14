module kofnga_types
  use kofnga_kinds, only : dp, i64
  implicit none
  private
  public :: kofnga_control, kofnga_result, kofnga_summary, objective_function

  abstract interface
    function objective_function(subset) result(value)
      import dp
      integer, intent(in) :: subset(:)
      real(dp) :: value
    end function objective_function
  end interface

  type :: kofnga_control
    integer :: popsize = 200
    integer :: keepbest = -1
    integer :: ngen = 500
    integer :: tourneysize = -1
    real(dp) :: mutprob = 0.01_dp
    real(dp) :: mutfrac = -1.0_dp
    integer(i64) :: seed = 5489_i64
    integer :: verbose = 0
  end type kofnga_control

  type :: kofnga_result
    integer, allocatable :: bestsol(:)
    real(dp) :: bestobj = huge(1.0_dp)
    integer, allocatable :: pop(:,:)
    real(dp), allocatable :: obj(:)
    integer, allocatable :: best_history(:,:)
    real(dp), allocatable :: obj_history(:)
    real(dp), allocatable :: avg_history(:)
  end type kofnga_result

  type :: kofnga_summary
    integer :: generations = 0
    integer :: unique_final = 0
    integer :: best_generation = 0
    integer, allocatable :: best_solution(:)
    real(dp) :: initial_average = 0.0_dp
    real(dp) :: initial_minimum = 0.0_dp
    real(dp) :: final_average = 0.0_dp
    real(dp) :: final_minimum = 0.0_dp
  end type kofnga_summary

end module kofnga_types
