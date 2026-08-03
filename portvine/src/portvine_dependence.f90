! SPDX-License-Identifier: GPL-3.0-only
module portvine_dependence
   use portvine_kinds, only : dp
   use portvine_types, only : asset_marginal_result, marginal_settings_type, &
      vine_settings_type, portvine_success, portvine_invalid_input, &
      portvine_vine_failure, vine_dvine, vine_cvine
   use portvine_marginals, only : get_marginal_point
   use portvine_ordering, only : greedy_dvine_order
   use rvinecopulib, only : dvine_model, cvine_model, make_dvine, make_cvine
   implicit none
   private
   public :: fit_vine_windows, vine_training_data

contains

   subroutine vine_training_data(marginal, marginal_window, first_index, n_train, data, status)
      type(asset_marginal_result), intent(in) :: marginal(:)
      integer, intent(in) :: marginal_window, first_index, n_train
      real(dp), intent(out) :: data(size(marginal),n_train)
      integer, intent(out), optional :: status
      real(dp) :: mu, sigma, u
      integer :: a, t, istat

      if (present(status)) status = portvine_success
      do a = 1, size(marginal)
         do t = 1, n_train
            call get_marginal_point(marginal(a),marginal_window,first_index+t-1, &
               mu,sigma,u,status=istat)
            if (istat /= portvine_success) then
               data = 0.5_dp
               if (present(status)) status = portvine_invalid_input
               return
            end if
            data(a,t) = u
         end do
      end do
   end subroutine vine_training_data

   subroutine fit_vine_windows(marginal, marginal_settings, vine_settings, nobs, &
      cond_indices, dvines, cvines, status)
      type(asset_marginal_result), intent(in) :: marginal(:)
      type(marginal_settings_type), intent(in) :: marginal_settings
      type(vine_settings_type), intent(in) :: vine_settings
      integer, intent(in) :: nobs
      integer, intent(in), optional :: cond_indices(:)
      type(dvine_model), allocatable, intent(out) :: dvines(:)
      type(cvine_model), allocatable, intent(out) :: cvines(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: data(:,:)
      integer, allocatable :: order(:)
      integer :: nv, v, mw, first_index, d, istat

      if (present(status)) status = portvine_success
      d = size(marginal)
      if (d < 2 .or. vine_settings%train_size < 3 .or. &
          vine_settings%refit_size < 1 .or. marginal_settings%refit_size < vine_settings%refit_size .or. &
          mod(marginal_settings%refit_size,vine_settings%refit_size) /= 0) then
         allocate(dvines(0),cvines(0))
         if (present(status)) status = portvine_invalid_input
         return
      end if
      nv = (nobs-marginal_settings%train_size+vine_settings%refit_size-1)/vine_settings%refit_size
      allocate(order(d),data(d,vine_settings%train_size))
      if (vine_settings%vine_type == vine_dvine) then
         allocate(dvines(nv),cvines(0))
      else if (vine_settings%vine_type == vine_cvine) then
         allocate(cvines(nv),dvines(0))
      else
         allocate(dvines(0),cvines(0))
         if (present(status)) status = portvine_invalid_input
         return
      end if
      do v = 1, nv
         mw = (v*vine_settings%refit_size+marginal_settings%refit_size-1)/ &
            marginal_settings%refit_size
         mw = max(1,min(size(marginal(1)%window),mw))
         first_index = marginal_settings%train_size+1+ &
            (v-1)*vine_settings%refit_size-vine_settings%train_size
         call vine_training_data(marginal,mw,first_index,vine_settings%train_size,data,istat)
         if (istat /= portvine_success) then
            if (present(status)) status = portvine_vine_failure
            return
         end if
         if (present(cond_indices)) then
            call greedy_dvine_order(data,order,cond_indices,vine_settings%cutoff_depth,istat)
         else
            call greedy_dvine_order(data,order,cutoff_depth=vine_settings%cutoff_depth,status=istat)
         end if
         if (istat /= 0) then
            if (present(status)) status = portvine_vine_failure
            return
         end if
         if (vine_settings%vine_type == vine_dvine) then
            dvines(v) = make_dvine(d,order)
            call fit_dvine_dispatch(dvines(v),data,vine_settings)
         else
            cvines(v) = make_cvine(d,order)
            call fit_cvine_dispatch(cvines(v),data,vine_settings)
         end if
      end do
   end subroutine fit_vine_windows

   subroutine fit_dvine_dispatch(model,data,settings)
      type(dvine_model), intent(inout) :: model
      real(dp), intent(in) :: data(:,:)
      type(vine_settings_type), intent(in) :: settings
      if (allocated(settings%families)) then
         if (size(settings%families) > 0) then
            call model%fit(data,settings%families,trim(settings%criterion),settings%allow_rotations)
         else
            call model%fit(data,criterion=trim(settings%criterion), &
               allow_rotations=settings%allow_rotations)
         end if
      else
         call model%fit(data,criterion=trim(settings%criterion), &
            allow_rotations=settings%allow_rotations)
      end if
   end subroutine fit_dvine_dispatch

   subroutine fit_cvine_dispatch(model,data,settings)
      type(cvine_model), intent(inout) :: model
      real(dp), intent(in) :: data(:,:)
      type(vine_settings_type), intent(in) :: settings
      if (allocated(settings%families)) then
         if (size(settings%families) > 0) then
            call model%fit(data,settings%families,trim(settings%criterion),settings%allow_rotations)
         else
            call model%fit(data,criterion=trim(settings%criterion), &
               allow_rotations=settings%allow_rotations)
         end if
      else
         call model%fit(data,criterion=trim(settings%criterion), &
            allow_rotations=settings%allow_rotations)
      end if
   end subroutine fit_cvine_dispatch

end module portvine_dependence
