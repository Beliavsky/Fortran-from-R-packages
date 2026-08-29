! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from R package nnet by W. N. Venables and B. D. Ripley.
module nnet_types
use r_compat, only: dp
implicit none
private
public :: nnet_model_t, multinom_model_t

type :: nnet_model_t
   integer :: n_inputs = 0
   integer :: n_hidden = 0
   integer :: n_outputs = 0
   integer :: n_units = 0
   integer :: first_hidden = 0
   integer :: first_output = 0
   integer :: ns_units = 0
   logical :: entropy = .false.
   logical :: softmax = .false.
   logical :: censored = .false.
   logical :: skip = .false.
   real(dp) :: value = 0.0_dp
   integer :: convergence = 0
   integer :: counts(2) = 0
   integer, allocatable :: nconn(:)   ! lower bound 0, stores 0-based cumulative counts
   integer, allocatable :: conn(:)    ! source unit index, bias = 0
   real(dp), allocatable :: wts(:)
   real(dp), allocatable :: decay(:)
   logical, allocatable :: mask(:)
   real(dp), allocatable :: fitted(:,:)
   real(dp), allocatable :: residuals(:,:)
   real(dp), allocatable :: hessian(:,:)
end type nnet_model_t

type :: multinom_model_t
   type(nnet_model_t) :: net
   integer :: n_classes = 0
   integer :: rank = 0
   integer :: edf = 0
   real(dp) :: deviance = 0.0_dp
   real(dp) :: aic = 0.0_dp
   real(dp), allocatable :: coefficients(:,:) ! (n_classes-1, ncol(x))
   real(dp), allocatable :: information(:,:)
   real(dp), allocatable :: covariance(:,:)
end type multinom_model_t

end module nnet_types
