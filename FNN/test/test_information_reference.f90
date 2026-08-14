program test_information_reference
  use fnn
  implicit none
  real(dp), parameter :: tol=3.0e-11_dp
  real(dp), target :: x(8,2), y(10,2)
  real(dp), allocatable :: h(:),ce(:),kl(:)
  real(dp) :: mi
  real(dp), parameter :: href(3)=[3.900026995310336_dp,3.2525335073166888_dp,3.135178627588788_dp]
  real(dp), parameter :: ceref(3)=[2.410016274351600_dp,2.3298730637363874_dp,2.659327258581673_dp]
  real(dp), parameter :: klref(3)=[-1.5029782807556378_dp,-0.9356280033772029_dp,-0.4888189288040165_dp]

  x(1,:)=[0.1_dp,0.4_dp]; x(2,:)=[0.7_dp,1.2_dp]; x(3,:)=[1.3_dp,0.2_dp]; x(4,:)=[1.8_dp,1.5_dp]
  x(5,:)=[2.4_dp,0.9_dp]; x(6,:)=[3.0_dp,2.2_dp]; x(7,:)=[3.7_dp,1.1_dp]; x(8,:)=[4.2_dp,2.8_dp]
  y(1,:)=[-0.2_dp,0.3_dp]; y(2,:)=[0.4_dp,1.0_dp]; y(3,:)=[1.0_dp,-0.1_dp]; y(4,:)=[1.5_dp,1.1_dp]
  y(5,:)=[2.0_dp,0.5_dp]; y(6,:)=[2.6_dp,1.8_dp]; y(7,:)=[3.1_dp,0.8_dp]; y(8,:)=[3.5_dp,2.4_dp]
  y(9,:)=[4.0_dp,1.6_dp]; y(10,:)=[4.6_dp,3.0_dp]

  h=entropy(x,3,"brute")
  ce=crossentropy(x,y,3,"brute")
  kl=kl_divergence(x,y,3,"brute")
  mi=mutinfo(x(:,1:1),x(:,2:2),2)
  call check(maxval(abs(h-href))<tol,"entropy reference")
  call check(maxval(abs(ce-ceref))<tol,"crossentropy reference")
  call check(maxval(abs(kl-klref))<tol,"KL reference")
  call check(abs(mi-0.007440476190476275_dp)<tol,"KSG MI reference")
  print *, "test_information_reference: PASS"
contains
  subroutine check(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg
    if(.not.ok) then
      print *, "FAIL: ",trim(msg)
      error stop 1
    end if
  end subroutine check
end program test_information_reference
