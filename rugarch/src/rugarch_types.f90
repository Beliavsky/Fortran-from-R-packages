! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_types
   use rugarch_kinds, only : dp
   use rugarch_distributions, only : dist_norm
   implicit none
   private

   integer, parameter, public :: model_sgarch    = 1
   integer, parameter, public :: model_gjrgarch  = 2
   integer, parameter, public :: model_egarch    = 3
   integer, parameter, public :: model_aparch    = 4
   integer, parameter, public :: model_igarch    = 5
   integer, parameter, public :: model_figarch   = 6
   integer, parameter, public :: model_csgarch   = 7
   integer, parameter, public :: model_realgarch = 8
   integer, parameter, public :: model_fgarch    = 9

   integer, parameter, public :: fgarch_garch    = 1
   integer, parameter, public :: fgarch_tgarch   = 2
   integer, parameter, public :: fgarch_avgarch  = 3
   integer, parameter, public :: fgarch_ngarch   = 4
   integer, parameter, public :: fgarch_nagarch  = 5
   integer, parameter, public :: fgarch_aparch   = 6
   integer, parameter, public :: fgarch_allgarch = 7
   integer, parameter, public :: fgarch_gjrgarch = 8

   ! Compatibility names used by the earlier fGarch translation.
   integer, parameter, public :: model_garch = model_sgarch

   type, public :: garch_spec
      integer :: model = model_sgarch
      integer :: cond_dist = dist_norm
      integer :: figarch_truncation = 500
      integer :: fgarch_submodel = fgarch_allgarch
      real(dp) :: mean = 0.0_dp
      real(dp) :: omega = 1.0e-6_dp
      real(dp) :: delta = 2.0_dp
      real(dp) :: skew = 1.0_dp
      real(dp) :: shape = 5.0_dp
      ! Generalized-hyperbolic distribution lambda.
      real(dp) :: lambda = 1.0_dp
      ! Hentschel fGARCH power lambda and asymmetry shifts.
      real(dp) :: fgarch_lambda = 2.0_dp
      real(dp) :: fgarch_fk = 1.0_dp
      real(dp) :: frac_d = 0.2_dp
      real(dp) :: rho = 0.95_dp
      real(dp) :: phi = 0.05_dp
      real(dp) :: xi = 0.0_dp
      real(dp) :: tau1 = 0.0_dp
      real(dp) :: tau2 = 0.0_dp
      real(dp) :: measurement_sd = 0.1_dp
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: eta1(:)
      real(dp), allocatable :: eta2(:)
   end type garch_spec

   type, public :: garch_fit_result
      type(garch_spec) :: spec
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 1
      character(len=160) :: message = 'not fitted'
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: measurement_residuals(:)
   end type garch_fit_result

   type, public :: distribution_fit_result
      real(dp) :: mean = 0.0_dp
      real(dp) :: sd = 1.0_dp
      real(dp) :: shape = 5.0_dp
      real(dp) :: skew = 1.0_dp
      real(dp) :: lambda = 1.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 1
   end type distribution_fit_result

   type, public :: forecast_result
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: sigma(:)
   end type forecast_result

   public :: make_garch_spec, model_name, configure_fgarch_submodel, fgarch_submodel_name

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
      allocate(spec%ar(0),spec%ma(0),spec%alpha(pp),spec%gamma(pp),spec%beta(qq), &
         spec%eta1(pp),spec%eta2(pp))
      if (pp > 0) then
         spec%alpha = 0.08_dp/real(pp,dp)
         spec%gamma = 0.0_dp
         spec%eta1 = 0.0_dp
         spec%eta2 = 0.0_dp
      end if
      if (qq > 0) spec%beta = 0.90_dp/real(qq,dp)
      if (spec%model == model_fgarch) call configure_fgarch_submodel(spec,fgarch_allgarch)
      if (spec%model == model_igarch .and. pp > 0 .and. qq > 0) then
         spec%alpha = 0.08_dp/real(pp,dp)
         spec%beta = 0.92_dp/real(qq,dp)
      end if
   end function make_garch_spec

   subroutine configure_fgarch_submodel(spec,submodel)
      type(garch_spec),intent(inout)::spec
      integer,intent(in)::submodel
      spec%fgarch_submodel=submodel
      select case(submodel)
      case(fgarch_garch)
         spec%fgarch_lambda=2.0_dp;spec%delta=2.0_dp;spec%fgarch_fk=0.0_dp
         if(allocated(spec%eta1))spec%eta1=0.0_dp
         if(allocated(spec%eta2))spec%eta2=0.0_dp
      case(fgarch_tgarch)
         spec%fgarch_lambda=1.0_dp;spec%delta=1.0_dp;spec%fgarch_fk=0.0_dp
         if(allocated(spec%eta2))spec%eta2=0.0_dp
      case(fgarch_avgarch)
         spec%fgarch_lambda=1.0_dp;spec%delta=1.0_dp;spec%fgarch_fk=0.0_dp
      case(fgarch_ngarch)
         spec%fgarch_lambda=2.0_dp;spec%delta=0.0_dp;spec%fgarch_fk=1.0_dp
         if(allocated(spec%eta1))spec%eta1=0.0_dp
         if(allocated(spec%eta2))spec%eta2=0.0_dp
      case(fgarch_nagarch)
         spec%fgarch_lambda=2.0_dp;spec%delta=2.0_dp;spec%fgarch_fk=0.0_dp
         if(allocated(spec%eta1))spec%eta1=0.0_dp
      case(fgarch_aparch)
         spec%fgarch_lambda=1.0_dp;spec%delta=0.0_dp;spec%fgarch_fk=1.0_dp
         if(allocated(spec%eta2))spec%eta2=0.0_dp
      case(fgarch_gjrgarch)
         spec%fgarch_lambda=2.0_dp;spec%delta=2.0_dp;spec%fgarch_fk=0.0_dp
         if(allocated(spec%eta2))spec%eta2=0.0_dp
      case default
         spec%fgarch_submodel=fgarch_allgarch
         spec%fgarch_lambda=2.0_dp;spec%delta=0.0_dp;spec%fgarch_fk=1.0_dp
      end select
   end subroutine configure_fgarch_submodel

   pure function fgarch_submodel_name(submodel) result(name)
      integer,intent(in)::submodel
      character(len=10)::name
      select case(submodel)
      case(fgarch_garch);name='GARCH'
      case(fgarch_tgarch);name='TGARCH'
      case(fgarch_avgarch);name='AVGARCH'
      case(fgarch_ngarch);name='NGARCH'
      case(fgarch_nagarch);name='NAGARCH'
      case(fgarch_aparch);name='APARCH'
      case(fgarch_gjrgarch);name='GJRGARCH'
      case default;name='ALLGARCH'
      end select
   end function fgarch_submodel_name

   pure function model_name(model) result(name)
      integer, intent(in) :: model
      character(len=12) :: name
      select case (model)
      case (model_sgarch);    name = 'sGARCH'
      case (model_gjrgarch);  name = 'gjrGARCH'
      case (model_egarch);    name = 'eGARCH'
      case (model_aparch);    name = 'apARCH'
      case (model_igarch);    name = 'iGARCH'
      case (model_figarch);   name = 'fiGARCH'
      case (model_csgarch);   name = 'csGARCH'
      case (model_realgarch); name = 'realGARCH'
      case (model_fgarch);    name = 'fGARCH'
      case default;           name = 'unknown'
      end select
   end function model_name

end module rugarch_types
