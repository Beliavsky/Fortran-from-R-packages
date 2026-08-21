program test_variance
  use gb2, only : dp,varscore_gb2,vepar_gb2,varscore_mixture,survey_score_variance
  use survey_types, only : survey_design_t
  use survey_design, only : make_design
  implicit none
  real(dp)::x(8),w(8),par(4),vsc(4,4),vc(4,4),u(6,2),vm(2,2),manual(2,2),ct(3,2),meanct(2)
  integer::cluster(6,1),i,j,k,fails
  type(survey_design_t)::des
  logical::ok
  x=[1.2_dp,1.8_dp,2.5_dp,3.1_dp,3.8_dp,4.6_dp,5.2_dp,6.5_dp]
  w=[1._dp,2._dp,1._dp,1.5_dp,.8_dp,1.2_dp,2.1_dp,.9_dp]
  par=[2.3_dp,4.2_dp,1.7_dp,3.4_dp]
  call varscore_gb2(x,par,vsc,w)
  call vepar_gb2(x,vsc,par,vc,w,ok=ok)
  fails=0
  if(.not.ok .or. maxval(abs(vc-transpose(vc)))>1e-10_dp .or. any([(vc(i,i)<0._dp,i=1,4)])) fails=fails+1
  u=reshape([1._dp,2._dp, 2._dp,1._dp, 3._dp,0._dp, 4._dp,1._dp, 5._dp,2._dp, 6._dp,1._dp],[6,2])
  cluster(:,1)=[1,1,2,2,3,3]
  call make_design([(1._dp,i=1,6)],cluster,des)
  call survey_score_variance(u,des,vm)
  ct(1,:)=u(1,:)+u(2,:)
  ct(2,:)=u(3,:)+u(4,:)
  ct(3,:)=u(5,:)+u(6,:)
  meanct=sum(ct,dim=1)/3._dp
  manual=0._dp
  do k=1,3
  do j=1,2
  do i=1,2
  manual(i,j)=manual(i,j)+1.5_dp*(ct(k,i)-meanct(i))*(ct(k,j)-meanct(j))
  end do
  end do
  end do
  if(maxval(abs(vm-manual))>1e-12_dp) then
  print *,'survey variance mismatch',vm,manual
  fails=fails+1
  end if
  call varscore_mixture(u,vm)
  if(any([(vm(i,i)<0._dp,i=1,2)])) fails=fails+1
  if(fails>0) error stop 1
  print '(a)','test_variance: PASS'
end program
