! SPDX-License-Identifier: GPL-2.0-or-later
program basic_wavelets
   use wavelets, only : dp, wavelet_transform_type, dwt_named, idwt
   implicit none

   type(wavelet_transform_type) :: wt
   real(dp) :: x(16, 1)
   real(dp), allocatable :: reconstructed(:, :)
   integer :: i
   integer :: ierr

   do i = 1, size(x, 1)
      x(i, 1) = real(i, dp)
   end do

   call dwt_named(x, wt, filter_name='haar', n_levels=3, ierr=ierr)
   if (ierr /= 0) error stop 'dwt failed'
   call idwt(wt, reconstructed, ierr)
   if (ierr /= 0) error stop 'idwt failed'

   write (*, '(a,i0)') 'levels: ', wt%level
   write (*, '(a,es12.4)') 'max reconstruction error: ', maxval(abs(reconstructed - x))
end program basic_wavelets
