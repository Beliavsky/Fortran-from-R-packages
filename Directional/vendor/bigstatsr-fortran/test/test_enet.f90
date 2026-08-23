program test_enet
    use bigstatsr
    implicit none
    type(fbm_real) :: x
    type(enet_path_result) :: eg,el
    real(dp) :: a(100,2),y(100),yl(100),lam(2),ctr(2),scl(2)
    integer :: i
    do i=1,100
        if(mod(i,2)==0) then
            a(i,1)=1.0_dp
        else
            a(i,1)=-1.0_dp
        end if
        a(i,2)=sin(0.3_dp*real(i,dp))
    end do
    ctr=[0.0_dp,sum(a(:,2))/100.0_dp]
    scl=[1.0_dp,sqrt(sum((a(:,2)-ctr(2))**2)/100.0_dp)]
    y=2.0_dp*a(:,1)
    yl=1.0_dp/(1.0_dp+exp(-(-0.2_dp+1.1_dp*a(:,1))))
    x=create_fbm('test_enet.bk',100,2)
    call fbm_from_array(x,a)
    lam=[0.25_dp,0.0_dp]
    eg=elastic_net_gaussian_path(x,y,lam,1.0_dp,ctr,scl,eps=1.0e-10_dp,maxiter=1000)
    call check(abs(eg%beta(1,2)-2.0_dp)<1.0e-8_dp,'gaussian lambda zero')
    call check(abs(eg%beta(1,1))<2.0_dp,'gaussian shrinkage')
    el=elastic_net_logistic_path(x,yl,lam,1.0_dp,ctr,scl,eps=1.0e-10_dp,maxiter=1000)
    call check(abs(el%beta(1,2)-1.1_dp)<2.0e-5_dp,'logistic lambda zero')
    call check(abs(el%intercept(2)+0.2_dp)<2.0e-5_dp,'logistic intercept')
    call execute_command_line('rm -f test_enet.bk')
    print *, 'test_enet: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok) then
            print *, 'FAIL: ',trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_enet
