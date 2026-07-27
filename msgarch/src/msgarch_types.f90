! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_types
   use msgarch_kinds, only : dp
   implicit none
   private
   type, public :: regime_spec
      character(len=12) :: model='sGARCH'
      character(len=8) :: distribution='norm'
      real(dp) :: omega=0.1_dp
      real(dp) :: alpha=0.1_dp
      real(dp) :: gamma=0.0_dp
      real(dp) :: beta=0.8_dp
      real(dp) :: shape=8.0_dp
      real(dp) :: skew=1.0_dp
   end type regime_spec

   type, public :: msgarch_spec
      integer :: k=0
      logical :: is_mixture=.false.
      type(regime_spec), allocatable :: regime(:)
      real(dp), allocatable :: transition(:,:)
   end type msgarch_spec

   type, public :: regime_state
      real(dp) :: h=1.0_dp
      real(dp) :: lnh=0.0_dp
      real(dp) :: fh=1.0_dp
   end type regime_state

   type, public :: filter_result
      real(dp) :: loglik=-huge(1.0_dp)
      real(dp), allocatable :: log_density(:,:)
      real(dp), allocatable :: variance(:,:)
      real(dp), allocatable :: predicted(:,:)
      real(dp), allocatable :: filtered(:,:)
      real(dp), allocatable :: smoothed(:,:)
      integer, allocatable :: viterbi(:)
      real(dp), allocatable :: next_probability(:)
   end type filter_result


   type, public :: posterior_state_result
      real(dp), allocatable :: predicted(:,:)
      real(dp), allocatable :: filtered(:,:)
      real(dp), allocatable :: smoothed(:,:)
   end type posterior_state_result

   type, public :: simulation_result
      real(dp), allocatable :: draw(:,:)
      integer, allocatable :: state(:,:)
      real(dp), allocatable :: conditional_sd(:,:,:)
   end type simulation_result

   type, public :: fit_result
      type(msgarch_spec) :: spec
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: hessian(:,:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: standard_error(:)
      real(dp) :: loglik=-huge(1.0_dp)
      real(dp) :: aic=huge(1.0_dp)
      real(dp) :: bic=huge(1.0_dp)
      integer :: iterations=0
      integer :: evaluations=0
      logical :: converged=.false.
   end type fit_result

   type, public :: mcmc_result
      real(dp), allocatable :: draws(:,:)
      real(dp), allocatable :: log_posterior(:)
      real(dp), allocatable :: posterior_mean(:)
      real(dp), allocatable :: posterior_sd(:)
      real(dp) :: acceptance_rate=0.0_dp
      real(dp) :: dic=huge(1.0_dp)
   end type mcmc_result

   type, public :: risk_result
      real(dp), allocatable :: var(:,:)
      real(dp), allocatable :: es(:,:)
      real(dp), allocatable :: simulations(:,:)
   end type risk_result

   type, public :: hmm_fit_result
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: transition(:,:)
      real(dp), allocatable :: probability(:)
      real(dp) :: loglik=-huge(1.0_dp)
      integer :: iterations=0
      logical :: converged=.false.
   end type hmm_fit_result
end module msgarch_types
