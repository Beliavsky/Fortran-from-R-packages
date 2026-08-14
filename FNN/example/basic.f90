program basic
  use fnn
  implicit none
  real(dp), target :: train(6,2),test(2,2)
  integer :: cl(6),i
  type(knn_result) :: z
  type(classification_result) :: pred

  train=reshape([0.0_dp,0.0_dp, 1.0_dp,0.0_dp, 0.0_dp,1.0_dp, &
                 4.0_dp,4.0_dp, 5.0_dp,4.0_dp, 4.0_dp,5.0_dp],[6,2],order=[2,1])
  test=reshape([0.2_dp,0.2_dp, 4.4_dp,4.2_dp],[2,2],order=[2,1])
  cl=[1,1,1,2,2,2]

  z=get_knnx(train,test,3,"kd_tree")
  pred=knn_classify(train,test,cl,3,"kd_tree")
  do i=1,size(test,1)
    print '(a,i0,a,*(i0,1x))',"query ",i," neighbor indices: ",z%index(i,:)
    print '(a,i0,a,f6.3)',"  predicted class ",pred%class(i),", probability = ",pred%probability(i)
  end do
end program basic
