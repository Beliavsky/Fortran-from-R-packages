! SPDX-License-Identifier: GPL-2.0-or-later
module esg
   use esg_models
   use esg_pricing
   use randompack, only: randompack_rng, randompack_rng_type
   implicit none
   public
contains
   pure subroutine esg_innovations(rng, innovations, status)
      type(randompack_rng_type), intent(inout) :: rng !! Caller-owned deterministic random stream; no global RNG is changed.
      real(dp), intent(out) :: innovations(:,:,:) !! Preallocated (scenarios,horizon,5) Gaussian innovations in risk-block order.
      integer, intent(out) :: status !! 0 success, 2 invalid shape, otherwise randompack normal-generation status.
      real(dp), allocatable :: draws(:)
      integer :: risk, n, h
      n = size(innovations, 1)
      h = size(innovations, 2)
      status = 2
      if (n < 1 .or. h < 1 .or. size(innovations, 3) /= 5) return
      allocate(draws(n * h))
      do risk = 1, 5
         call rng%normal(draws, info=status)
         if (status /= 0) return
         innovations(:,:,risk) = reshape(draws, [n, h])
      end do
   end subroutine
end module
