! Computational translation of qrng 0.0-11 generalized Halton generator.
! Upstream qrng license: GPL-2 | GPL-3.
module qrng_ghalton_mod
use r_mod, only: dp, runif1
use qrng_halton_data_mod, only: ghalton_max_dim, ghalton_primes, ghalton_perm_tn2
implicit none
private
public :: ghalton

contains

function ghalton(n, d, method, shift_coeff) result(u)
integer, intent(in) :: n, d
character(len=*), intent(in), optional :: method
integer, intent(in), optional :: shift_coeff(:,:)
real(dp), allocatable :: u(:,:)
integer, allocatable :: shift(:,:)
integer :: base, coeff(32), f
integer :: i, j, k, tmp
real(dp) :: value
logical :: generalized
character(len=:), allocatable :: meth

if (n < 1) error stop "ghalton: n must be at least 1"
if (d < 1 .or. d > ghalton_max_dim) error stop "ghalton: d must be in 1,...,360"

meth = "generalized"
if (present(method)) meth = trim(adjustl(method))
select case (meth)
case ("generalized")
   generalized = .true.
case ("halton")
   generalized = .false.
case default
   error stop "ghalton: method must be 'generalized' or 'halton'"
end select

allocate(u(n,d), shift(d,32))
if (present(shift_coeff)) then
   if (size(shift_coeff,1) /= d .or. size(shift_coeff,2) /= 32) &
      error stop "ghalton: shift_coeff must have shape (d,32)"
   shift = shift_coeff
   do j = 1, d
      if (any(shift(j,:) < 0) .or. any(shift(j,:) >= ghalton_primes(j))) &
         error stop "ghalton: shift coefficients must be base digits"
   end do
else
   do j = 1, d
      base = ghalton_primes(j)
      do k = 1, 32
         shift(j,k) = int(real(base,dp) * runif1())
      end do
   end do
end if

! Upstream qrng randomizes generalized and plain Halton by a 32-digit shift.
do j = 1, d
   base = ghalton_primes(j)
   value = 0.0_dp
   do k = 32, 1, -1
      value = (value + real(shift(j,k),dp)) / real(base,dp)
   end do
   u(1,j) = value
end do

do i = 1, n - 1
   do j = 1, d
      base = ghalton_primes(j)
      coeff = 0
      tmp = i
      k = 1
      do while (tmp > 0 .and. k <= 32)
         coeff(k) = modulo(tmp, base)
         tmp = tmp / base
         k = k + 1
      end do
      if (generalized) then
         f = ghalton_perm_tn2(j)
      else
         f = 1
      end if
      value = 0.0_dp
      do k = 32, 1, -1
         value = value + real(modulo(f * coeff(k) + shift(j,k), base), dp)
         value = value / real(base, dp)
      end do
      u(i+1,j) = value
   end do
end do
end function ghalton

end module qrng_ghalton_mod
