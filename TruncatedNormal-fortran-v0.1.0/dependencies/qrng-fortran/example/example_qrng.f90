program example_qrng
use qrng, only: dp, sobol
implicit none
real(dp), allocatable :: u(:,:)
integer :: i

u = sobol(8,3)
do i = 1, size(u,1)
   print '(3f10.6)', u(i,:)
end do
end program example_qrng
