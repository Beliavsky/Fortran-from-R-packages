! Computational test functions translated from qrng 0.0-11 R/test_functions.R.
! Upstream qrng license: GPL-2 | GPL-3.
module qrng_test_functions_mod
use r_compat, only: dp, sort
implicit none
private
public :: sum_of_squares, sobol_g, exceedance_indicator
public :: exceedance_individual_given_sum, exceedance_sum_given_sum

contains

pure function sum_of_squares(u) result(value)
real(dp), intent(in) :: u(:,:)
real(dp), allocatable :: value(:)
if (size(u,2) < 1) error stop "sum_of_squares: u must have at least one column"
allocate(value(size(u,1)))
value = 3.0_dp * sum(u*u, dim=2) / real(size(u,2),dp)
end function sum_of_squares

pure function sobol_g(u, alpha) result(value)
! Independent-copula kernel of qrng::sobol_g().  For a non-independent copula,
! apply its inverse Rosenblatt transform before calling this function.
real(dp), intent(in) :: u(:,:)
real(dp), intent(in), optional :: alpha(:)
real(dp), allocatable :: value(:), a(:)
integer :: i, j, d
real(dp) :: term
d = size(u,2)
if (d < 1) error stop "sobol_g: u must have at least one column"
allocate(a(d), value(size(u,1)))
if (present(alpha)) then
   if (size(alpha) /= d) error stop "sobol_g: alpha must have length size(u,2)"
   a = alpha
else
   do j = 1, d
      a(j) = real(j,dp)
   end do
end if
value = 1.0_dp
do i = 1, size(u,1)
   do j = 1, d
      term = (abs(4.0_dp*u(i,j)-2.0_dp) + a(j)) / (1.0_dp + a(j))
      value(i) = value(i) * term
   end do
end do
end function sobol_g

function exceedance_indicator(x, q, p) result(indicator)
real(dp), intent(in) :: x(:,:)
real(dp), intent(in), optional :: q(:), p(:)
logical, allocatable :: indicator(:)
real(dp), allocatable :: threshold(:), probs(:)
integer :: i, j, d
d = size(x,2)
if (d < 1) error stop "exceedance_indicator: x must have at least one column"
allocate(threshold(d), indicator(size(x,1)))
if (present(q)) then
   if (size(q) == 1) then
      threshold = q(1)
   else if (size(q) == d) then
      threshold = q
   else
      error stop "exceedance_indicator: q must have length 1 or d"
   end if
else
   allocate(probs(d))
   probs = 0.99_dp
   if (present(p)) then
      if (size(p) == 1) then
         probs = p(1)
      else if (size(p) == d) then
         probs = p
      else
         error stop "exceedance_indicator: p must have length 1 or d"
      end if
   end if
   if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp)) &
      error stop "exceedance_indicator: p must lie in [0,1]"
   do j = 1, d
      threshold(j) = empirical_quantile_type1(x(:,j), probs(j))
   end do
end if
indicator = .true.
do i = 1, size(x,1)
   do j = 1, d
      if (.not. (x(i,j) > threshold(j))) then
         indicator(i) = .false.
         exit
      end if
   end do
end do
end function exceedance_indicator

function exceedance_individual_given_sum(x, q, p) result(out)
real(dp), intent(in) :: x(:,:)
real(dp), intent(in), optional :: q, p
real(dp), allocatable :: out(:,:)
real(dp), allocatable :: s(:)
logical, allocatable :: keep(:)
real(dp) :: threshold, prob
integer :: i, k
allocate(s(size(x,1)))
s = sum(x,dim=2)
prob = 0.99_dp
if (present(p)) prob = p
if (prob < 0.0_dp .or. prob > 1.0_dp) error stop "exceedance: p must lie in [0,1]"
if (present(q)) then
   threshold = q
else
   threshold = empirical_quantile_type1(s,prob)
end if
keep = s > threshold
allocate(out(count(keep),size(x,2)))
k = 0
do i = 1, size(x,1)
   if (keep(i)) then
      k = k + 1
      out(k,:) = x(i,:)
   end if
end do
end function exceedance_individual_given_sum

function exceedance_sum_given_sum(x, q, p) result(out)
real(dp), intent(in) :: x(:,:)
real(dp), intent(in), optional :: q, p
real(dp), allocatable :: out(:), s(:)
real(dp) :: threshold, prob
allocate(s(size(x,1)))
s = sum(x,dim=2)
prob = 0.99_dp
if (present(p)) prob = p
if (prob < 0.0_dp .or. prob > 1.0_dp) error stop "exceedance: p must lie in [0,1]"
if (present(q)) then
   threshold = q
else
   threshold = empirical_quantile_type1(s,prob)
end if
out = pack(s,s > threshold)
end function exceedance_sum_given_sum

pure real(dp) function empirical_quantile_type1(x,p) result(q)
! R quantile(..., type=1): inverse empirical distribution function.
real(dp), intent(in) :: x(:), p
real(dp), allocatable :: xs(:)
integer :: idx, n
if (size(x) < 1) error stop "empirical_quantile_type1: empty input"
if (p < 0.0_dp .or. p > 1.0_dp) error stop "empirical_quantile_type1: p must lie in [0,1]"
xs = sort(x)
n = size(xs)
idx = max(1, min(n, ceiling(real(n,dp)*p)))
q = xs(idx)
end function empirical_quantile_type1

end module qrng_test_functions_mod
