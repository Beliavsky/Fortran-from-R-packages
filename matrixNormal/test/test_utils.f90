program test_utils
  use matrixNormal
  use test_support
  implicit none
  real(dp),allocatable::a(:,:),b(:,:),v(:),h(:),k(:,:)
  logical::ok
  a=reshape([1.0_dp,2.0_dp,2.0_dp,4.0_dp],[2,2])
  call assert_close(tr(a),5.0_dp,1e-14_dp,'trace')
  v=vec(a)
  call assert_true(all(abs(v-[1.0_dp,2.0_dp,2.0_dp,4.0_dp])<1e-14_dp),'vec column order')
  h=vech(a,ok=ok)
  call assert_true(ok,'vech symmetric status')
  call assert_true(all(abs(h-[1.0_dp,2.0_dp,4.0_dp])<1e-14_dp),'vech lower triangle')
  call assert_true(is_square_matrix(a),'square matrix')
  call assert_true(is_symmetric_matrix(a),'symmetric matrix')
  call assert_true(is_positive_semidefinite(a),'positive semidefinite matrix')
  call assert_true(.not.is_positive_definite(a),'singular matrix is not positive definite')
  b=reshape([2.0_dp,-1.0_dp,-1.0_dp,2.0_dp],[2,2])
  call assert_true(is_positive_definite(b),'positive definite matrix')
  call assert_true(all(abs(identity_matrix(3)-reshape([1.0_dp,0.0_dp,0.0_dp, &
    0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp],[3,3]))<1e-14_dp), &
    'identity_matrix')
  call assert_true(all(abs(ones_matrix(2,3)-1.0_dp)<1e-14_dp),'ones_matrix')
  k=kronecker_product(reshape([1.0_dp,3.0_dp,2.0_dp,4.0_dp],[2,2]),reshape([5.0_dp,7.0_dp,6.0_dp,8.0_dp],[2,2]))
  call assert_close(k(1,1),5.0_dp,1e-14_dp,'kron first')
  call assert_close(k(4,4),32.0_dp,1e-14_dp,'kron last')
  write(*,'(a)') 'test_utils: PASS'
end program test_utils
