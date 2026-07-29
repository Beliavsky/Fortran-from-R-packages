! SPDX-License-Identifier: MIT
program test_matrix
  use bekks
  implicit none
  real(dp) :: a(3,3),v(6),x(6,2),ylag_expected(4,5)
  real(dp), allocatable :: ylag(:,:)
  integer :: i
  a=reshape([1.0_dp,2.0_dp,3.0_dp,2.0_dp,4.0_dp,5.0_dp,3.0_dp,5.0_dp,6.0_dp],[3,3])
  v=vech_lower(a)
  if(maxval(abs(v-[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]))>1.0e-14_dp)error stop 'vech'
  if(maxval(abs(matmul(elimination_mat(3),reshape(a,[9]))-v))>1.0e-14_dp)error stop 'elimination'
  if(maxval(abs(matmul(duplication_mat(3),v)-reshape(a,[9])))>1.0e-14_dp)error stop 'duplication'
  if(maxval(abs(matmul(commutation_mat(3),reshape(a,[9]))-reshape(transpose(a),[9])))>1.0e-14_dp)error stop 'commutation'
  do i=1,6;x(i,:)=[real(i,dp),10.0_dp+real(i,dp)];end do
  ylag=y_lag_cr(x,2)
  ylag_expected=reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
    2.0_dp,3.0_dp,4.0_dp,5.0_dp,12.0_dp,13.0_dp,14.0_dp,15.0_dp, &
    1.0_dp,2.0_dp,3.0_dp,4.0_dp,11.0_dp,12.0_dp,13.0_dp,14.0_dp],[4,5])
  if(maxval(abs(ylag-ylag_expected))>1.0e-14_dp)error stop 'y lag'
  print '(a)','test_matrix: PASS'
end program test_matrix
