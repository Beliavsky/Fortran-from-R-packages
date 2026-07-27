! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_types
   use gogarch_kinds, only : dp
   implicit none
   private

   type, public :: univariate_spec
      character(len=12) :: model = 'garch'
      character(len=12) :: distribution = 'norm'
      integer :: p = 1
      integer :: o = 0
      integer :: q = 1
      logical :: include_mean = .true.
      logical :: fit_delta = .true.
      logical :: fit_shape = .true.
      logical :: fit_skew = .true.
      real(dp) :: delta = 2.0_dp
      real(dp) :: shape = 8.0_dp
      real(dp) :: skew = 1.0_dp
   end type univariate_spec

   type, public :: garch11_fit
      character(len=12) :: model = 'garch'
      character(len=12) :: distribution = 'norm'
      integer :: p = 1
      integer :: o = 0
      integer :: q = 1
      real(dp) :: mean = 0.0_dp
      real(dp) :: omega = 0.0_dp
      real(dp) :: alpha = 0.0_dp
      real(dp) :: beta = 0.0_dp
      real(dp) :: delta = 2.0_dp
      real(dp) :: shape = 8.0_dp
      real(dp) :: skew = 1.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 1
      real(dp), allocatable :: arch(:)
      real(dp), allocatable :: leverage(:)
      real(dp), allocatable :: garch(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: power_scale(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: standardized(:)
   end type garch11_fit

   type, public :: gogarch_fit
      character(len=16) :: method = ''
      integer :: n = 0
      integer :: m = 0
      integer :: status = 1
      integer :: optimizer_iterations = 0
      real(dp) :: log_likelihood = -huge(1.0_dp)
      type(univariate_spec) :: factor_spec
      real(dp), allocatable :: data(:,:)
      real(dp), allocatable :: sample_covariance(:,:)
      real(dp), allocatable :: covariance_sqrt(:,:)
      real(dp), allocatable :: covariance_invsqrt(:,:)
      real(dp), allocatable :: eigenvectors(:,:)
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: rotation(:,:)
      real(dp), allocatable :: mixing(:,:)
      real(dp), allocatable :: factors(:,:)
      type(garch11_fit), allocatable :: factor_models(:)
      real(dp), allocatable :: factor_variance(:,:)
      real(dp), allocatable :: covariance(:,:,:)
      real(dp), allocatable :: mm_weights(:)
      real(dp), allocatable :: objective_parameters(:)
   end type gogarch_fit

end module gogarch_types
