! SPDX-License-Identifier: GPL-2.0-or-later
program transform_example
   use moments, only : dp, all_moments, raw2central, central2raw, all_cumulants
   implicit none

   real(dp) :: x(6)
   real(dp), allocatable :: raw(:), central(:), recovered(:), cumulants(:)
   integer :: k

   x = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 5.0_dp]
   raw = all_moments(x, 5)
   central = raw2central(raw)
   recovered = central2raw(central, raw(2))
   cumulants = all_cumulants(raw)

   write(*, '(a,es12.4)') 'round-trip maximum error: ', maxval(abs(recovered - raw))
   do k = 0, 5
      write(*, '(a,i0,a,f14.6)') 'cumulant ', k, ': ', cumulants(k + 1)
   end do
end program transform_example
