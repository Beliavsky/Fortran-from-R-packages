! FNN-fortran: modern Fortran translation of computational code from FNN 1.1.4.1.
! Modified/translated 2026 by the FNN-fortran contributors.
! SPDX-License-Identifier: GPL-2.0-or-later
! See UPSTREAM.md and upstream/FNN-1.1.4.1 for original authorship and notices.
module fnn_types
  use fnn_kinds, only : dp
  implicit none
  private
  public :: knn_result, classification_result, regression_result, ownn_result

  type :: knn_result
    integer, allocatable :: index(:,:)
    real(dp), allocatable :: distance(:,:)
  end type knn_result

  type :: classification_result
    integer, allocatable :: class(:)
    real(dp), allocatable :: probability(:)
    integer, allocatable :: nn_index(:,:)
    real(dp), allocatable :: nn_distance(:,:)
  end type classification_result

  type :: regression_result
    real(dp), allocatable :: prediction(:)
    real(dp), allocatable :: residuals(:)
    real(dp) :: press = 0.0_dp
    real(dp) :: r2_predict = 0.0_dp
    logical :: cross_validated = .false.
  end type regression_result

  type :: ownn_result
    integer :: k = 0
    integer, allocatable :: knn_class(:)
    integer, allocatable :: ownn_class(:)
    integer, allocatable :: bnn_class(:)
    real(dp), allocatable :: knn_probability(:)
    real(dp), allocatable :: ownn_probability(:)
    real(dp), allocatable :: bnn_probability(:)
    real(dp) :: accuracy(3) = 0.0_dp
    logical :: has_accuracy = .false.
  end type ownn_result
end module fnn_types
