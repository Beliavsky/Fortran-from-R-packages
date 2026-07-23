! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

module tseries_fft
   use tseries_kinds, only : dp
   implicit none
   private

   public :: discrete_fourier_transform
   public :: inverse_discrete_fourier_transform

contains

   subroutine discrete_fourier_transform(x,z)
      real(dp), intent(in) :: x(:)
      complex(dp), intent(out) :: z(:)
      integer :: n,j,k
      real(dp), parameter :: two_pi=2.0_dp*acos(-1.0_dp)
      complex(dp) :: phase
      n=size(x)
      if(size(z)/=n) return
      do k=1,n
         z(k)=(0.0_dp,0.0_dp)
         do j=1,n
            phase=cmplx(0.0_dp,-two_pi*real((j-1)*(k-1),dp)/real(n,dp),dp)
            z(k)=z(k)+x(j)*exp(phase)
         end do
      end do
   end subroutine discrete_fourier_transform

   subroutine inverse_discrete_fourier_transform(z,x)
      complex(dp), intent(in) :: z(:)
      real(dp), intent(out) :: x(:)
      integer :: n,j,k
      real(dp), parameter :: two_pi=2.0_dp*acos(-1.0_dp)
      complex(dp) :: value,phase
      n=size(z)
      if(size(x)/=n) return
      do j=1,n
         value=(0.0_dp,0.0_dp)
         do k=1,n
            phase=cmplx(0.0_dp,two_pi*real((j-1)*(k-1),dp)/real(n,dp),dp)
            value=value+z(k)*exp(phase)
         end do
         x(j)=real(value,dp)/real(n,dp)
      end do
   end subroutine inverse_discrete_fourier_transform

end module tseries_fft
