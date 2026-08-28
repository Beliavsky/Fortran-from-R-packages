! Computational translation of qrng 0.0-11 Sobol generator.
! Upstream qrng license: GPL-2 | GPL-3.
module qrng_sobol_mod
use, intrinsic :: iso_fortran_env, only: int64
use r_mod, only: dp, runif1, set_seed_int
use qrng_sobol_data_mod, only: sobol_max_dim, sobol_max_degree, sobol_poly, sobol_minit
implicit none
private
public :: sobol

contains

function sobol(n, d, randomize, skip, seed, digital_shift) result(u)
integer, intent(in) :: n, d
logical, intent(in), optional :: randomize
integer, intent(in), optional :: skip, seed
real(dp), intent(in), optional :: digital_shift(:)
real(dp), allocatable :: u(:,:)
integer(int64), allocatable :: v(:,:), lastpoint(:)
real(dp), allocatable :: rvector(:)
integer(int64) :: total, maxn, gray, count, point, randint
integer(int64), parameter :: rmaxint = 4503599627370496_int64 ! 2^52
integer(int64) :: newv
integer :: numcols, i, j, k, degree, temp, column, sk
real(dp) :: recipd, urand, scaled
logical :: do_randomize

if (n < 1) error stop "sobol: n must be at least 1"
if (d < 1 .or. d > sobol_max_dim) error stop "sobol: d must be in 1,...,16510"
if (int(n,int64) > 2147483647_int64) error stop "sobol: n must be <= 2^31-1"
sk = 0
if (present(skip)) sk = skip
if (sk < 0) error stop "sobol: skip must be nonnegative"
total = int(n,int64) + int(sk,int64)
if (total > 2147483648_int64) &
   error stop "sobol: n + skip must be <= 2^31 for the 32-column upstream table"

do_randomize = .false.
if (present(randomize)) do_randomize = randomize
if (present(seed)) call set_seed_int(seed)
if (present(digital_shift)) then
   if (size(digital_shift) /= d) error stop "sobol: digital_shift must have length d"
   if (any(digital_shift < 0.0_dp) .or. any(digital_shift >= 1.0_dp)) &
      error stop "sobol: digital_shift values must lie in [0,1)"
end if

numcols = 0
maxn = 1_int64
do while (maxn < total)
   maxn = maxn * 2_int64
   numcols = numcols + 1
end do
if (numcols > 31) error stop "sobol: requested index requires more than 31 active columns"

allocate(u(n,d), lastpoint(d))
allocate(v(d,max(1,numcols)))
v = 0_int64
lastpoint = 0_int64

if (numcols > 0) then
   v(1,1:numcols) = 1_int64
   do i = 2, d
      degree = sobol_max_degree
      do while (degree > 0)
         if (btest(int(sobol_poly(i),int64), degree)) exit
         degree = degree - 1
      end do
      do j = 1, min(degree,numcols)
         v(i,j) = int(sobol_minit(j,i-1), int64)
      end do
      do j = degree + 1, numcols
         newv = v(i,j-degree)
         temp = 1
         do k = degree - 1, 0, -1
            if (btest(int(sobol_poly(i),int64), k)) then
               newv = ieor(newv, shiftl(v(i,j-(degree-k)), temp))
            end if
            temp = temp + 1
         end do
         v(i,j) = newv
      end do
   end do

   temp = 1
   do j = numcols - 1, 1, -1
      do i = 1, d
         v(i,j) = shiftl(v(i,j), temp)
      end do
      temp = temp + 1
   end do
end if

recipd = 1.0_dp / real(maxn,dp)
if (do_randomize) then
   allocate(rvector(d))
   if (present(digital_shift)) then
      rvector = digital_shift
   else
      do i = 1, d
         do
            urand = runif1()
            scaled = urand * real(maxn,dp)
            if (real(floor(scaled,kind=int64),dp) /= scaled) exit
         end do
         rvector(i) = urand
      end do
   end if
end if

if (sk > 0 .and. numcols > 0) then
   gray = ieor(int(sk,int64), shiftr(int(sk,int64),1))
   do j = 0, numcols - 1
      if (btest(gray,j)) then
         do i = 1, d
            lastpoint(i) = ieor(lastpoint(i), v(i,j+1))
         end do
      end if
   end do
end if

call store_point(1)

do count = int(sk,int64), int(sk,int64) + int(n,int64) - 2_int64
   column = 0
   do while (btest(count,column))
      column = column + 1
   end do
   do i = 1, d
      lastpoint(i) = ieor(lastpoint(i), v(i,column+1))
   end do
   call store_point(int(count - int(sk,int64) + 2_int64))
end do

contains

subroutine store_point(row)
integer, intent(in) :: row
integer :: ii
if (do_randomize) then
   do ii = 1, d
      point = shiftl(lastpoint(ii), 52 - numcols)
      randint = floor(rvector(ii) * real(rmaxint,dp), kind=int64)
      point = ieor(point, randint)
      u(row,ii) = real(point,dp) / real(rmaxint,dp)
   end do
else
   do ii = 1, d
      u(row,ii) = real(lastpoint(ii),dp) * recipd
   end do
end if
end subroutine store_point

end function sobol

end module qrng_sobol_mod
