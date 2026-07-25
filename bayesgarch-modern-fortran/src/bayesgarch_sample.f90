! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
module bayesgarch_sample
   use bayesgarch_kinds, only : dp
   use bayesgarch_sampler, only : bayesgarch_result
   implicit none
   private

   public :: form_posterior_sample
   public :: posterior_mean
   public :: posterior_sd

contains

   subroutine form_posterior_sample(result, burn, thin, sample)
      type(bayesgarch_result), intent(in) :: result
      integer, intent(in) :: burn
      integer, intent(in) :: thin
      real(dp), allocatable, intent(out) :: sample(:, :)
      integer :: n_chains
      integer :: n_iter
      integer :: n_keep
      integer :: chain
      integer :: iteration
      integer :: row

      if (.not. allocated(result%draws)) error stop "form_posterior_sample: result has no draws"
      n_iter = size(result%draws, 1)
      n_chains = size(result%draws, 3)
      if (burn < 0 .or. burn >= n_iter) error stop "form_posterior_sample: invalid burn"
      if (thin < 1) error stop "form_posterior_sample: thin must be positive"
      n_keep = (n_iter - burn - 1) / thin + 1
      allocate(sample(n_keep * n_chains, 4))
      row = 0
      do chain = 1, n_chains
         do iteration = burn + 1, n_iter, thin
            row = row + 1
            sample(row, :) = result%draws(iteration, :, chain)
         end do
      end do
   end subroutine form_posterior_sample

   pure function posterior_mean(sample) result(value)
      real(dp), intent(in) :: sample(:, :)
      real(dp) :: value(size(sample, 2))

      value = sum(sample, dim=1) / real(size(sample, 1), dp)
   end function posterior_mean

   pure function posterior_sd(sample) result(value)
      real(dp), intent(in) :: sample(:, :)
      real(dp) :: value(size(sample, 2))
      real(dp) :: mean(size(sample, 2))
      integer :: i

      if (size(sample, 1) < 2) then
         value = 0.0_dp
         return
      end if
      mean = posterior_mean(sample)
      value = 0.0_dp
      do i = 1, size(sample, 1)
         value = value + (sample(i, :) - mean)**2
      end do
      value = sqrt(value / real(size(sample, 1) - 1, dp))
   end function posterior_sd

end module bayesgarch_sample
