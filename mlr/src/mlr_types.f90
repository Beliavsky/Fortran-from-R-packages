module mlr_types
  use mlr_kinds, only : dp
  implicit none
  private
  public :: index_set, resample_plan, scaler_model, linear_model, logistic_model
  public :: kmeans_model, knn_model, metric_summary, tune_result, feature_select_result
  type :: index_set
    integer, allocatable :: idx(:)
  end type index_set

  type :: resample_plan
    type(index_set), allocatable :: train(:), test(:)
    integer, allocatable :: group(:)
  end type resample_plan

  type :: scaler_model
    real(dp), allocatable :: center(:), scale(:)
  end type scaler_model

  type :: linear_model
    real(dp), allocatable :: coef(:)
    real(dp) :: sigma = 0.0_dp
    logical :: fitted = .false.
  end type linear_model

  type :: logistic_model
    real(dp), allocatable :: coef(:)
    integer :: iterations = 0
    logical :: converged = .false.
  end type logistic_model

  type :: kmeans_model
    real(dp), allocatable :: centers(:,:)
    integer :: iterations = 0
    real(dp) :: withinss = 0.0_dp
    logical :: converged = .false.
  end type kmeans_model

  type :: knn_model
    real(dp), allocatable :: x(:,:)
    real(dp), allocatable :: y(:)
    integer, allocatable :: class_y(:)
    integer :: k = 1
    logical :: classification = .false.
  end type knn_model

  type :: metric_summary
    real(dp) :: mean = 0.0_dp
    real(dp) :: sd = 0.0_dp
    real(dp), allocatable :: values(:)
  end type metric_summary

  type :: tune_result
    real(dp), allocatable :: par(:)
    real(dp) :: objective = huge(1.0_dp)
    integer :: evaluations = 0
  end type tune_result

  type :: feature_select_result
    logical, allocatable :: selected(:)
    real(dp) :: objective = huge(1.0_dp)
    integer :: evaluations = 0
  end type feature_select_result
end module mlr_types
