program test_pspline
  use survival
  implicit none
  real(dp) :: x(5)=[0._dp,0.25_dp,0.5_dp,0.75_dp,1._dp]
  real(dp),allocatable :: b(:,:),p(:,:),knots(:)
  integer :: ierr
  call pspline_basis(x,8,3,[0._dp,1._dp],.false.,b,p,knots,ierr)
  if(ierr/=0) error stop 'pspline status'
  if(size(b,1)/=5 .or. size(b,2)/=10) error stop 'pspline shape'
  if(maxval(abs(p-transpose(p)))>1e-12_dp) error stop 'pspline symmetry'
  print *, 'test_pspline PASS'
end program
