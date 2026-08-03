! SPDX-License-Identifier: GPL-3.0-only
module portvine_types
   use portvine_kinds, only : dp
   use rugarch, only : garch_spec, garch_fit_result, make_garch_spec, &
      model_sgarch, dist_sstd
   use rvinecopulib, only : dvine_model, cvine_model
   implicit none
   private

   integer, parameter, public :: portvine_success = 0
   integer, parameter, public :: portvine_invalid_input = 1
   integer, parameter, public :: portvine_fit_failure = 2
   integer, parameter, public :: portvine_vine_failure = 3

   integer, parameter, public :: risk_var = 1
   integer, parameter, public :: risk_es_mean = 2
   integer, parameter, public :: risk_es_median = 3
   integer, parameter, public :: risk_es_mc = 4

   integer, parameter, public :: vine_dvine = 1
   integer, parameter, public :: vine_cvine = 2
   integer, parameter, public :: vine_rvine_approx = vine_cvine

   type, public :: marginal_settings_type
      integer :: train_size = 0
      integer :: refit_size = 0
      type(garch_spec), allocatable :: spec(:)
      integer :: max_iterations = 1800
   end type marginal_settings_type

   type, public :: vine_settings_type
      integer :: train_size = 0
      integer :: refit_size = 0
      integer :: vine_type = vine_dvine
      integer, allocatable :: families(:)
      character(len=8) :: criterion = 'aic'
      logical :: allow_rotations = .true.
      integer :: cutoff_depth = huge(1)
   end type vine_settings_type

   type, public :: arma_fit_result
      real(dp) :: intercept = 0.0_dp
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: css = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 1
   end type arma_fit_result

   type, public :: marginal_window_result
      integer :: start_index = 0
      integer :: forecast_start = 0
      integer :: forecast_end = 0
      type(garch_fit_result) :: garch_fit
      type(garch_spec) :: combined_spec
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: standardized(:)
      real(dp), allocatable :: uniform(:)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: sigma(:)
      integer :: status = 1
   end type marginal_window_result

   type, public :: asset_marginal_result
      type(marginal_window_result), allocatable :: window(:)
      integer :: status = 1
   end type asset_marginal_result

   type, public :: portvine_roll_result
      real(dp), allocatable :: overall(:,:,:)
      real(dp), allocatable :: conditional(:,:,:,:)
      real(dp), allocatable :: realized(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: condition_level(:)
      real(dp), allocatable :: conditional_value(:,:,:)
      integer, allocatable :: row_index(:)
      integer, allocatable :: vine_window(:)
      integer, allocatable :: measure(:)
      real(dp), allocatable :: weights(:,:)
      type(asset_marginal_result), allocatable :: marginal(:)
      type(dvine_model), allocatable :: dvine(:)
      type(cvine_model), allocatable :: cvine(:)
      integer :: vine_type = vine_dvine
      integer :: status = portvine_success
      character(len=240) :: message = 'ok'
   end type portvine_roll_result

   public :: make_portvine_spec, make_marginal_settings, make_vine_settings
   public :: risk_measure_name

contains

   function make_portvine_spec(ar_order, ma_order, arch_order, garch_order, &
      model, cond_dist) result(spec)
      integer, intent(in), optional :: ar_order, ma_order, arch_order, garch_order
      integer, intent(in), optional :: model, cond_dist
      type(garch_spec) :: spec
      integer :: pa, qm, p, q, imod, idist

      pa = 1
      qm = 1
      p = 1
      q = 1
      imod = model_sgarch
      idist = dist_sstd
      if (present(ar_order)) pa = max(0, ar_order)
      if (present(ma_order)) qm = max(0, ma_order)
      if (present(arch_order)) p = max(0, arch_order)
      if (present(garch_order)) q = max(0, garch_order)
      if (present(model)) imod = model
      if (present(cond_dist)) idist = cond_dist
      spec = make_garch_spec(p, q, imod, idist)
      if (allocated(spec%ar)) deallocate(spec%ar)
      if (allocated(spec%ma)) deallocate(spec%ma)
      allocate(spec%ar(pa), spec%ma(qm))
      if (pa > 0) spec%ar = 0.0_dp
      if (qm > 0) spec%ma = 0.0_dp
   end function make_portvine_spec

   function make_marginal_settings(train_size, refit_size, n_assets, default_spec) result(settings)
      integer, intent(in) :: train_size, refit_size, n_assets
      type(garch_spec), intent(in), optional :: default_spec
      type(marginal_settings_type) :: settings
      type(garch_spec) :: spec
      integer :: j

      settings%train_size = train_size
      settings%refit_size = refit_size
      allocate(settings%spec(max(0,n_assets)))
      spec = make_portvine_spec()
      if (present(default_spec)) spec = default_spec
      do j = 1, n_assets
         settings%spec(j) = spec
      end do
   end function make_marginal_settings

   function make_vine_settings(train_size, refit_size, vine_type, families) result(settings)
      integer, intent(in) :: train_size, refit_size
      integer, intent(in), optional :: vine_type, families(:)
      type(vine_settings_type) :: settings

      settings%train_size = train_size
      settings%refit_size = refit_size
      if (present(vine_type)) settings%vine_type = vine_type
      if (present(families)) then
         allocate(settings%families(size(families)))
         settings%families = families
      else
         allocate(settings%families(0))
      end if
   end function make_vine_settings

   pure function risk_measure_name(measure) result(name)
      integer, intent(in) :: measure
      character(len=12) :: name
      select case (measure)
      case (risk_var)
         name = 'VaR'
      case (risk_es_mean)
         name = 'ES_mean'
      case (risk_es_median)
         name = 'ES_median'
      case (risk_es_mc)
         name = 'ES_mc'
      case default
         name = 'unknown'
      end select
   end function risk_measure_name

end module portvine_types
