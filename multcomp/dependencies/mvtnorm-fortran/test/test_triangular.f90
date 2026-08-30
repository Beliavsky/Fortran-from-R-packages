! SPDX-License-Identifier: GPL-2.0-only
program test_triangular
  use mvtnorm
  use test_support
  implicit none
  real(dp) :: cov(3,3),eye(3,3)
  real(dp),allocatable :: l(:,:),linv(:,:),back(:,:),cor(:,:),pre(:,:),pc(:,:),v(:),u(:,:),score(:)
  logical :: ok
  character(len=256) :: message

  cov=reshape([2.0_dp,0.6_dp,-0.2_dp, 0.6_dp,1.5_dp,0.3_dp, -0.2_dp,0.3_dp,1.0_dp],[3,3])
  l=cov2chol(cov,ok,message); call assert_true(ok,'cov2chol')
  back=chol2cov(l); call assert_matrix_close(back,cov,2.0e-14_dp,'chol2cov')
  linv=chol2invchol(l,ok); call assert_true(ok,'chol2invchol')
  eye=matmul(l,linv); call assert_matrix_close(eye,identity_matrix(3),2.0e-14_dp,'triangular inverse')
  back=invchol2cov(linv); call assert_matrix_close(back,cov,3.0e-14_dp,'invchol2cov')
  pre=invchol2pre(linv); back=matmul(pre,cov); call assert_matrix_close(back,identity_matrix(3),5.0e-14_dp,'precision')
  cor=chol2cor(l); call assert_close(cor(1,1),1.0_dp,1.0e-15_dp,'cor diagonal')
  call assert_close(cor(1,2),0.6_dp/sqrt(3.0_dp),2.0e-14_dp,'correlation')
  pc=chol2pc(l); call assert_true(maxval(abs(pc-transpose(pc)))<=0.0_dp,'partial correlation symmetry')
  v=pack_lower(l,.true.,.false.); u=unpack_lower(v,3,.true.,.false.)
  call assert_matrix_close(u,l,0.0_dp,'pack/unpack')
  score=deperma_score(l,[1,2,3],spread(1.0_dp,1,size(v)),1.0e-6_dp)
  call assert_close(maxval(abs(score-1.0_dp)),0.0_dp,2.0e-8_dp,'identity depermutation score')
  print '(a)', 'test_triangular: PASS'
end program test_triangular
