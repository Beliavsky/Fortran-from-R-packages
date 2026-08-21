program test_simplex
    use boot_kinds, only : dp
    use boot_simplex
    implicit none
    real(dp)::c(2),a(3,2),b(3)
    type(simplex_result)::res
    c=[3.0_dp,2.0_dp]
    a=reshape([1.0_dp,1.0_dp, 1.0_dp,0.0_dp, 0.0_dp,1.0_dp],[3,2],order=[2,1])
    ! Explicit rows after reshape are not portable to read mentally; overwrite.
    a(1,:)=[1.0_dp,1.0_dp]
    a(2,:)=[1.0_dp,0.0_dp]
    a(3,:)=[0.0_dp,1.0_dp]
    b=[4.0_dp,2.0_dp,3.0_dp]
    call simplex_solve(c,res,a_le=a,b_le=b,maximize=.true.)
    if(res%status/=1)error stop 1
    if(abs(res%value-10.0_dp)>1.0e-8_dp)error stop 2
    if(maxval(abs(res%x-[2.0_dp,2.0_dp]))>1.0e-8_dp)error stop 3
    print '(a)', 'test_simplex: PASS'
end program test_simplex
