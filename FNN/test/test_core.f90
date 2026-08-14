program test_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fnn
  implicit none
  real(dp), parameter :: tol=2.0e-12_dp
  real(dp), target :: x(6,2), q(2,2), yreg(6), xc(4,2), xdup(4,1)
  integer :: cl(6)
  type(knn_result) :: a,b,c,cr
  type(classification_result) :: cls,cv
  type(regression_result) :: reg
  type(ownn_result) :: ow
  real(dp), allocatable :: h(:),ce(:),kl(:),mie(:)
  real(dp) :: mi

  x=reshape([0.0_dp,0.0_dp, 1.0_dp,0.1_dp, 2.0_dp,0.0_dp, &
             0.0_dp,2.0_dp, 2.0_dp,2.1_dp, 4.0_dp,4.0_dp],[6,2],order=[2,1])
  q=reshape([0.2_dp,0.1_dp, 3.0_dp,3.0_dp],[2,2],order=[2,1])
  a=get_knn(x,3,"brute")
  b=get_knn(x,3,"kd_tree")
  c=get_knn(x,3,"cover_tree")
  call assert_true(all(a%index==b%index),"kd self indices")
  call assert_true(all(a%index==c%index),"cover self indices")
  call assert_true(maxval(abs(a%distance-b%distance))<tol,"kd self distances")
  call assert_true(maxval(abs(a%distance-c%distance))<tol,"cover self distances")

  a=get_knnx(x,q,3,"brute"); b=get_knnx(x,q,3,"kd_tree"); c=get_knnx(x,q,3,"cover_tree")
  call assert_true(all(a%index==b%index),"kd query indices")
  call assert_true(all(a%index==c%index),"cover query indices")
  call assert_true(maxval(abs(a%distance-b%distance))<tol,"kd query distances")
  call assert_true(maxval(abs(a%distance-c%distance))<tol,"cover query distances")


  xdup(:,1)=[0.0_dp,0.0_dp,1.0_dp,2.0_dp]
  a=get_knn(xdup,2,"kd_tree")
  call assert_true(all(a%index(:,1)/=[1,2,3,4]),"duplicate self exclusion")
  c=get_knn(xdup,2,"cover_tree")
  call assert_true(all(a%index==c%index),"duplicate cover agreement")

  xc=reshape([1.0_dp,0.0_dp, 0.0_dp,1.0_dp, -1.0_dp,0.0_dp, 0.0_dp,-1.0_dp],[4,2],order=[2,1])
  cr=get_knn(xc,2,"CR")
  call assert_true(cr%index(1,1)==2,"CR tie order")
  call assert_close(cr%distance(1,1),1.0_dp,tol,"CR distance")
  call assert_close(cr%distance(1,2),1.0_dp,tol,"CR distance 2")

  h=entropy(x,2,"brute")
  ce=crossentropy(x+0.35_dp,x,2,"kd_tree")
  kl=kl_divergence(x,x+0.35_dp,2,"kd_tree")
  call assert_true(all(ieee_is_finite(h)),"entropy finite")
  call assert_true(all(ieee_is_finite(ce)),"crossentropy finite")
  call assert_true(all(ieee_is_finite(kl)),"KL finite")
  call assert_true(size(h)==2 .and. size(ce)==2 .and. size(kl)==2,"information vector sizes")

  mi=mutinfo(x(:,1:1),x(:,2:2),2)
  mie=mutual_information_entropy(x(:,1:1),x(:,2:2),2,"brute")
  call assert_true(ieee_is_finite(mi),"mutinfo finite")
  call assert_true(size(mie)==2,"entropy MI size")

  cl=[1,1,1,2,2,2]
  cls=knn_classify(x,q,cl,3,"kd_tree")
  call assert_true(cls%class(1)==1 .and. cls%class(2)==2,"classification")
  call assert_true(all(cls%probability>0.0_dp),"classification probability")
  cv=knn_cv(x,cl,1,"cover_tree")
  call assert_true(size(cv%class)==6,"classification CV")

  yreg=[0.0_dp,1.0_dp,2.0_dp,2.0_dp,4.0_dp,8.0_dp]
  reg=knn_reg(x,yreg,2,test=q,algorithm="kd_tree")
  call assert_true(size(reg%prediction)==2,"regression prediction")
  reg=knn_reg(x,yreg,2,algorithm="brute")
  call assert_true(reg%cross_validated,"regression CV flag")
  call assert_true(size(reg%residuals)==6,"regression residuals")

  ow=ownn(x,q,cl,k=3,algorithm="kd_tree",testcl=[1,2])
  call assert_true(ow%k==3,"ownn k")
  call assert_true(size(ow%ownn_class)==2 .and. size(ow%bnn_class)==2,"ownn outputs")
  call assert_true(ow%has_accuracy,"ownn accuracy")
  ow=ownn(x,q,cl,algorithm="brute",seed=42)
  call assert_true(ow%k>=1 .and. ow%k<=size(x,1),"ownn automatic k")

  print *, "test_core: PASS"
contains
  subroutine assert_true(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg
    if(.not.ok) then
      print *, "FAIL: ",trim(msg)
      error stop 1
    end if
  end subroutine assert_true
  subroutine assert_close(a,b,t,msg)
    real(dp), intent(in) :: a,b,t
    character(len=*), intent(in) :: msg
    call assert_true(abs(a-b)<=t,msg)
  end subroutine assert_close
end program test_core
