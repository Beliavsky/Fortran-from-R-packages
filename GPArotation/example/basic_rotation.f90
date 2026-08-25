program basic_rotation
  use gpa_kinds, only: dp
  use gpa_rotation
  use gpa_criteria, only: criterion_options
  implicit none
  real(dp) :: a(6,2)
  type(rotation_result) :: r
  type(rotation_options) :: opt
  type(criterion_options) :: co
  integer :: i
  a(1,:)=[0.80_dp,0.25_dp]; a(2,:)=[0.75_dp,0.30_dp]; a(3,:)=[0.70_dp,0.20_dp]
  a(4,:)=[0.25_dp,0.80_dp]; a(5,:)=[0.30_dp,0.72_dp]; a(6,:)=[0.20_dp,0.68_dp]
  opt%eps=1.0e-7_dp
  call gpforth(a,'varimax',r,co,opt)
  print '(a,l1)', 'converged: ',r%converged
  print '(a,f12.6)', 'objective: ',r%objective
  print '(a)', 'rotated loadings:'
  do i=1,size(a,1); print '(2f11.6)',r%loadings(i,:); end do
end program basic_rotation
