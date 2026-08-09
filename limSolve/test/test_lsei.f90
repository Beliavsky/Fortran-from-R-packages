program test_lsei
    use limsolve
    implicit none
    type(solve_result) :: r
    real(dp) :: a(2,2),b(2),e0(0,2),f0(0),g(2,2),h(2)
    real(dp) :: e(1,2),f(1)
    real(dp) :: au(4,3),bu(4),eu(1,3),fu(1),gu(2,3),hu(2)

    a=0.0_dp; a(1,1)=1.0_dp; a(2,2)=1.0_dp
    b=[2.0_dp,-1.0_dp]
    g=0.0_dp; g(1,1)=1.0_dp; g(2,2)=1.0_dp; h=0.0_dp
    call lsei(a,b,e0,f0,g,h,r)
    call check(r%succeeded(),'lsei nonnegative status')
    call check(maxval(abs(r%x-[2.0_dp,0.0_dp])) < 2.0e-5_dp,'lsei nonnegative solution')

    e(1,:)=[1.0_dp,1.0_dp]; f=[1.0_dp]
    b=0.0_dp
    call lsei(a,b,e,f,g,h,r,fulloutput=.true.)
    call check(r%succeeded(),'lsei equality status')
    call check(maxval(abs(r%x-[0.5_dp,0.5_dp])) < 2.0e-5_dp,'lsei equality solution')
    call check(r%rank_eq==1,'lsei rank')
    call check(allocated(r%covariance),'lsei covariance')

    call ldei(e,f,g,h,r)
    call check(r%succeeded(),'ldei status')
    call check(maxval(abs(r%x-[0.5_dp,0.5_dp])) < 2.0e-6_dp,'ldei solution')


    au=reshape([3.0_dp,1.0_dp,2.0_dp,0.0_dp,2.0_dp,0.0_dp,0.0_dp,1.0_dp, &
        1.0_dp,0.0_dp,2.0_dp,0.0_dp],[4,3])
    bu=[2.0_dp,1.0_dp,8.0_dp,3.0_dp]
    eu(1,:)=[0.0_dp,1.0_dp,0.0_dp]; fu=[3.0_dp]
    gu=reshape([-1.0_dp,1.0_dp,2.0_dp,0.0_dp,0.0_dp,-1.0_dp],[2,3])
    hu=[-3.0_dp,2.0_dp]
    call lsei(au,bu,eu,fu,gu,hu,r)
    call check(r%succeeded(),'upstream lsei example status')
    call check(maxval(abs(r%x-[1.2424242424_dp,3.0_dp,-0.7575757576_dp])) < 2.0e-7_dp, &
        'upstream lsei example solution')
    call check(abs(r%solution_norm-98.0606060606_dp) < 2.0e-7_dp, &
        'upstream lsei example objective')

    print *, 'PASS test_lsei'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(len=*),intent(in)::msg
        if(.not.ok) then; print *, 'FAIL: ',trim(msg); error stop 1; end if
    end subroutine check
end program test_lsei
