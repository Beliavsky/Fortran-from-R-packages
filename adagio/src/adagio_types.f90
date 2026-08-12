! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of computational code from R package adagio 0.9.2.
module adagio_types
  use adagio_kinds, only : dp
  implicit none
  private
  public :: opt_result, ea_result, de_result, cmaes_result
  public :: assignment_result, subset_result, knapsack_result
  public :: mknapsack_result, setcover_result, change_result
  public :: binpack_result, maxsub_result, maxsub2d_result, maxempty_result
  public :: hamiltonian_result, count_result

  type :: opt_result
     real(dp), allocatable :: x(:)
     real(dp) :: f = huge(1.0_dp)
     integer :: evaluations = 0
     integer :: iterations = 0
     integer :: restarts = 0
     integer :: convergence = 0
  end type

  type :: ea_result
     real(dp), allocatable :: par(:)
     real(dp) :: val = huge(1.0_dp)
     integer :: fun_calls = 0
     integer :: actual_fun_calls = 0
     real(dp) :: rel_scl = huge(1.0_dp)
     real(dp) :: rel_tol = huge(1.0_dp)
  end type

  type :: de_result
     real(dp), allocatable :: xmin(:)
     real(dp) :: fmin = huge(1.0_dp)
     integer :: nfeval = 0
     integer :: actual_nfeval = 0
  end type

  type :: cmaes_result
     real(dp), allocatable :: xmin(:)
     real(dp) :: fmin = huge(1.0_dp)
     integer :: evaluations = 0
  end type

  type :: assignment_result
     integer, allocatable :: perm(:)
     real(dp) :: value = 0.0_dp
     integer :: err = 0
  end type

  type :: subset_result
     integer, allocatable :: inds(:)
     integer :: val = 0
     logical :: found = .false.
  end type

  type :: knapsack_result
     integer, allocatable :: indices(:)
     integer :: capacity = 0
     real(dp) :: profit = 0.0_dp
  end type

  type :: mknapsack_result
     integer, allocatable :: ksack(:)
     real(dp) :: value = 0.0_dp
     integer :: bs = 0
  end type

  type :: setcover_result
     integer, allocatable :: sets(:)
     real(dp) :: objective = huge(1.0_dp)
     logical :: feasible = .false.
  end type

  type :: change_result
     integer :: count = 0
     integer, allocatable :: solution(:)
     logical :: feasible = .false.
  end type

  type :: binpack_result
     integer :: nbins = 0
     integer, allocatable :: xbins(:)
     real(dp), allocatable :: sbins(:)
     real(dp) :: filled = 0.0_dp
  end type

  type :: maxsub_result
     real(dp) :: sum = 0.0_dp
     integer :: first = 0
     integer :: last = 0
  end type

  type :: maxsub2d_result
     real(dp) :: sum = 0.0_dp
     integer :: inds(4) = 0
     real(dp), allocatable :: submat(:,:)
  end type

  type :: maxempty_result
     real(dp) :: area = 0.0_dp
     real(dp) :: rect(4) = 0.0_dp
  end type

  type :: hamiltonian_result
     integer, allocatable :: path(:)
     logical :: found = .false.
  end type

  type :: count_result
     real(dp), allocatable :: values(:)
     integer, allocatable :: counts(:)
  end type
end module adagio_types
