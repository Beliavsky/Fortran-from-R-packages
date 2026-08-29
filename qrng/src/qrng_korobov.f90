! Computational translation of qrng 0.0-11 Korobov generator.
! Upstream qrng license: GPL-2 | GPL-3.
module qrng_korobov_mod
use, intrinsic :: iso_fortran_env, only: int64
use r_compat, only: dp, runif1
implicit none
private
public :: korobov

interface korobov
   module procedure korobov_scalar_generator
   module procedure korobov_vector_generator
end interface korobov

contains

function korobov_scalar_generator(n, d, generator, randomize) result(u)
integer, intent(in) :: n, d, generator
logical, intent(in), optional :: randomize
real(dp), allocatable :: u(:,:)
integer :: g(1)
g(1) = generator
u = korobov_vector_generator(n, d, g, randomize)
end function korobov_scalar_generator

function korobov_vector_generator(n, d, generator, randomize) result(u)
integer, intent(in) :: n, d
integer, intent(in) :: generator(:)
logical, intent(in), optional :: randomize
real(dp), allocatable :: u(:,:)
integer(int64), allocatable :: g(:)
real(dp), allocatable :: aux(:)
real(dp) :: shift
integer :: i, j
logical :: do_randomize

if (n < 2) error stop "korobov: n must be at least 2"
if (d < 1) error stop "korobov: d must be at least 1"
if (int(n, int64) > 2147483647_int64) error stop "korobov: n must be <= 2^31-1"
if (int(d, int64) > 2147483647_int64) error stop "korobov: d must be <= 2^31-1"
if (size(generator) /= 1 .and. size(generator) /= d) &
   error stop "korobov: generator must have length 1 or d"
if (any(generator < 1) .or. any(generator > n - 1)) &
   error stop "korobov: generator entries must be in 1,...,n-1"

do_randomize = .false.
if (present(randomize)) do_randomize = randomize

allocate(g(d), aux(d), u(n,d))
if (size(generator) == 1) then
   do j = 1, d
      g(j) = modular_power(int(generator(1), int64), j - 1, int(n, int64))
   end do
else
   g = int(generator, int64)
end if

aux = real(g, dp) / real(n, dp)
u(1,:) = 0.0_dp
do i = 2, n
   do j = 1, d
      u(i,j) = u(i-1,j) + aux(j)
      if (u(i,j) > 1.0_dp) u(i,j) = u(i,j) - 1.0_dp
   end do
end do

if (do_randomize) then
   do j = 1, d
      shift = runif1()
      do i = 1, n
         u(i,j) = u(i,j) + shift
         if (u(i,j) > 1.0_dp) u(i,j) = u(i,j) - 1.0_dp
      end do
   end do
end if
end function korobov_vector_generator

pure integer(int64) function modular_power(base, exponent, modulus) result(value)
integer(int64), intent(in) :: base, modulus
integer, intent(in) :: exponent
integer(int64) :: b
integer :: e
value = 1_int64
b = modulo(base, modulus)
e = exponent
do while (e > 0)
   if (btest(e, 0)) value = modulo(value * b, modulus)
   e = shiftr(e, 1)
   if (e > 0) b = modulo(b * b, modulus)
end do
value = modulo(value, modulus)
end function modular_power

end module qrng_korobov_mod
