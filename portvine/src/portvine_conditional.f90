! SPDX-License-Identifier: GPL-3.0-only
module portvine_conditional
   use portvine_kinds, only : dp
   use portvine_types, only : portvine_success, portvine_invalid_input
   use rvinecopulib, only : dvine_model
   implicit none
   private
   public :: conditional_dvine_sample

contains

   subroutine conditional_dvine_sample(model, n_samples, conditioning, quantile_mode, sample, status)
      type(dvine_model), intent(in) :: model
      integer, intent(in) :: n_samples
      real(dp), intent(in) :: conditioning(:)
      logical, intent(in) :: quantile_mode
      real(dp), intent(out) :: sample(model%dimension,n_samples)
      integer, intent(out), optional :: status
      real(dp), allocatable :: z(:,:)
      real(dp) :: u
      integer :: d, i, j, nc

      d = model%dimension
      nc = size(conditioning)
      if (present(status)) status = portvine_success
      if (d < 2 .or. n_samples < 1 .or. nc < 1 .or. nc > 2 .or. nc > d .or. &
          any(conditioning <= 0.0_dp) .or. any(conditioning >= 1.0_dp)) then
         sample = 0.0_dp
         if (present(status)) status = portvine_invalid_input
         return
      end if
      allocate(z(d,n_samples))
      do j = 1, n_samples
         do i = 1, d
            call random_number(u)
            z(i,j) = max(1.0e-10_dp,min(1.0_dp-1.0e-10_dp,u))
         end do
      end do
      z(model%order(1),:) = conditioning(1)
      if (nc == 2) then
         if (quantile_mode) then
            ! The second supplied value is a conditional quantile level.
            z(model%order(2),:) = conditioning(2)
         else
            ! Convert the raw second copula value to the Rosenblatt scale.
            z(model%order(2),:) = model%pair(1,2)%hfunc1(conditioning(1),conditioning(2))
         end if
      end if
      call model%inverse_rosenblatt(z,sample)
   end subroutine conditional_dvine_sample

end module portvine_conditional
