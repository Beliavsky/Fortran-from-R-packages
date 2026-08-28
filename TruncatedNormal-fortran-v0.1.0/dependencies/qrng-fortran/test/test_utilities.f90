program test_utilities
use qrng, only: dp, to_array_matrix, to_array_3d, sum_of_squares, sobol_g, &
   exceedance_indicator, exceedance_individual_given_sum, exceedance_sum_given_sum
implicit none
real(dp) :: x(2,6), z(5,2), u(2,2)
real(dp), allocatable :: y(:,:), a(:,:,:), v(:), rows(:,:), sums(:)
logical, allocatable :: ind(:)
integer :: i

x(1,:) = [1,2,3,4,5,6]
x(2,:) = [7,8,9,10,11,12]
y = to_array_matrix(x,3)
if (any(shape(y) /= [6,2])) error stop "to_array_matrix shape"
if (any(abs(y(:,1)-[1,3,5,7,9,11]) > 0.0_dp)) error stop "to_array_matrix col1"
if (any(abs(y(:,2)-[2,4,6,8,10,12]) > 0.0_dp)) error stop "to_array_matrix col2"
a = to_array_3d(x,3)
if (any(shape(a) /= [2,3,2])) error stop "to_array_3d shape"
if (any(abs(a(:,:,1)-reshape([1.0_dp,7.0_dp,2.0_dp,8.0_dp,3.0_dp,9.0_dp],[2,3])) > 0.0_dp)) &
   error stop "to_array_3d first slab"

u = reshape([0.0_dp,0.5_dp,1.0_dp,0.25_dp],[2,2])
v = sum_of_squares(u)
if (any(abs(v - 1.5_dp*sum(u*u,dim=2)) > 1.0e-14_dp)) error stop "sum_of_squares"
v = sobol_g(u, [1.0_dp,2.0_dp])
if (any(v <= 0.0_dp)) error stop "sobol_g"

do i=1,5
   z(i,:) = real(i,dp)
end do
ind = exceedance_indicator(z,p=[0.6_dp])
if (any(ind .neqv. [.false.,.false.,.false.,.true.,.true.])) error stop "exceedance_indicator"
rows = exceedance_individual_given_sum(z,p=0.6_dp)
if (any(shape(rows) /= [2,2])) error stop "exceedance rows shape"
if (any(abs(rows(:,1)-[4.0_dp,5.0_dp]) > 0.0_dp)) error stop "exceedance rows"
sums = exceedance_sum_given_sum(z,p=0.6_dp)
if (any(abs(sums-[8.0_dp,10.0_dp]) > 0.0_dp)) error stop "exceedance sums"

print '(a)', 'test_utilities: PASS'
end program test_utilities
