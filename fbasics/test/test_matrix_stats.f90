! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
program test_matrix_stats
  use fbasics
  use test_support
  implicit none
  real(dp),allocatable::h(:,:),p(:,:),ainv(:,:),kronm(:,:),pd(:,:),lags(:,:),z(:,:)
  real(dp)::a(2,2),x(5),mat(3,4),lm1,lm2,t3,t4,rm(3)
  type(basic_stats_result)::bs
  integer::info
  a=reshape([4.0_dp,2.0_dp,2.0_dp,3.0_dp],[2,2])
  call matrix_inverse(a,ainv,info)
  call assert_true(info==0,'matrix inverse info')
  call assert_close(ainv(1,1),0.375_dp,1e-12_dp,'matrix inverse')
  h=hilbert_matrix(3);call assert_close(h(3,2),0.25_dp,1e-14_dp,'Hilbert matrix')
  p=pascal_matrix(4);call assert_close(p(4,4),20.0_dp,1e-14_dp,'Pascal matrix')
  call assert_true(matrix_rank(a)==2,'matrix rank')
  call assert_close(matrix_norm(a,1),6.0_dp,1e-12_dp,'one norm')
  call assert_close(matrix_trace(a),7.0_dp,1e-14_dp,'trace')
  kronm=kronecker_product(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2]),reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]))
  call assert_close(kronm(3,3),4.0_dp,1e-14_dp,'Kronecker product')
  call assert_true(is_positive_definite(a),'positive definite check')
  call make_positive_definite(reshape([1.0_dp,2.0_dp,2.0_dp,1.0_dp],[2,2]),pd,info=info)
  call assert_true(info==0,'positive definite repair info')
  call assert_true(is_positive_definite(pd),'positive definite repair')
  x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
  lags=lag_matrix(x,[-1,0,2],.true.)
  call assert_true(size(lags,1)==2.and.size(lags,2)==3,'trimmed lag dimensions')
  call assert_close(lags(1,1),4.0_dp,1e-14_dp,'lead value')
  z=polynomial_distributed_lags(x,1,2,.true.)
  call assert_true(size(z,1)==3.and.size(z,2)==2,'PDL dimensions')
  bs=basic_stats(x)
  call assert_close(bs%mean,3.0_dp,1e-14_dp,'basic mean')
  call assert_close(bs%variance,2.5_dp,1e-14_dp,'basic variance')
  mat=reshape([(real(info,dp),info=1,12)],[3,4])
  rm=row_means(mat)
  call assert_close(rm(2),6.5_dp,1e-14_dp,'row mean')
  call sample_lmoments(x,lm1,lm2,t3,t4)
  call assert_close(lm1,3.0_dp,1e-14_dp,'L mean')
  call assert_true(lm2>0.0_dp,'L scale')
  call set_lcg_seed(12345_8)
  call assert_close(runif_lcg(),0.09661652850760917_dp,1e-14_dp,'LCG first value')
  call assert_close(heaviside(0.0_dp),0.5_dp,1e-14_dp,'Heaviside zero')
  write(*,'(a)')'Matrix, lag, RNG, and statistics tests passed.'
end program
