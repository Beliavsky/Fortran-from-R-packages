! Part of the experimental modern Fortran translation of fGarch 4052.93.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original fGarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

module fgarch_types
   use fgarch_kinds, only : dp
   use fgarch_distributions, only : dist_norm
   implicit none
   private

   integer, parameter, public :: model_garch = 1
   integer, parameter, public :: model_aparch = 2
   integer, parameter, public :: model_egarch = 3

   type, public :: garch_spec
      integer :: model = model_garch
      integer :: cond_dist = dist_norm
      real(dp) :: mean = 0.0_dp
      real(dp) :: omega = 1.0e-6_dp
      real(dp) :: delta = 2.0_dp
      real(dp) :: skew = 1.0_dp
      real(dp) :: shape = 5.0_dp
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: beta(:)
   end type garch_spec

   type, public :: garch_fit_result
      type(garch_spec) :: spec
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 1
      character(len=128) :: message = 'not fitted'
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: sigma(:)
   end type garch_fit_result

   type, public :: distribution_fit_result
      real(dp) :: mean = 0.0_dp
      real(dp) :: sd = 1.0_dp
      real(dp) :: shape = 5.0_dp
      real(dp) :: skew = 1.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 1
   end type distribution_fit_result

   public :: make_garch_spec

contains

   function make_garch_spec(p, q, model, cond_dist) result(spec)
      integer, intent(in), optional :: p, q, model, cond_dist
      type(garch_spec) :: spec
      integer :: pp, qq

      pp = 1
      qq = 1
      if (present(p)) pp = max(0,p)
      if (present(q)) qq = max(0,q)
      if (present(model)) spec%model = model
      if (present(cond_dist)) spec%cond_dist = cond_dist
      allocate(spec%ar(0),spec%ma(0),spec%alpha(pp),spec%gamma(pp),spec%beta(qq))
      if (pp > 0) then
         spec%alpha = 0.1_dp/real(pp,dp)
         spec%gamma = 0.0_dp
      end if
      if (qq > 0) spec%beta = 0.8_dp/real(qq,dp)
      if (spec%model == model_garch) spec%delta = 2.0_dp
   end function make_garch_spec

end module fgarch_types
