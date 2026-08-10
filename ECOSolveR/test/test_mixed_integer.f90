program test_mixed_integer
    use ecos
    implicit none
    type(ecos_dims)::d
    type(ecos_result)::r
    real(dp)::c2(2),g2(2,2),h2(2),g4(4,2),h4(4)
    real(dp)::c6(6),g6(3,6),h6(3)
    integer::bv1(1),iv2(2),bv6(6)

    d%l=2; c2=[-1.1_dp,-1.0_dp]
    g2=reshape([2.0_dp,3.0_dp,1.0_dp,4.0_dp],[2,2]); h2=[4.0_dp,12.0_dp]; bv1=[1]
    call ecos_csolve(c2,g2,h2,d,r,bool_vars=bv1)
    call check(r%exitflag==ECOS_OPTIMAL,'bool status')
    call check(maxval(abs(r%x-[1.0_dp,2.0_dp]))<1.0e-6_dp,'bool solution')

    d%l=4; c2=[-1.0_dp,-1.1_dp]
    g4=reshape([2.0_dp,3.0_dp,-1.0_dp,0.0_dp,1.0_dp,4.0_dp,0.0_dp,-1.0_dp],[4,2])
    h4=[4.0_dp,12.0_dp,0.0_dp,0.0_dp]; iv2=[1,2]
    call ecos_csolve(c2,g4,h4,d,r,int_vars=iv2)
    call check(r%exitflag==ECOS_OPTIMAL,'integer status')
    call check(maxval(abs(r%x-[0.0_dp,3.0_dp]))<1.0e-6_dp,'integer solution')

    d%l=3
    c6=[3.0_dp,5.0_dp,6.0_dp,9.0_dp,10.0_dp,10.0_dp]
    g6=reshape([2.0_dp,5.0_dp,-5.0_dp,-6.0_dp,3.0_dp,1.0_dp,3.0_dp,-1.0_dp,-4.0_dp, &
        -4.0_dp,-3.0_dp,2.0_dp,-1.0_dp,2.0_dp,-2.0_dp,2.0_dp,-1.0_dp,1.0_dp],[3,6])
    h6=[-2.0_dp,2.0_dp,-3.0_dp]; bv6=[1,2,3,4,5,6]
    call ecos_csolve(c6,g6,h6,d,r,bool_vars=bv6)
    call check(r%exitflag==ECOS_OPTIMAL,'six-binary status')
    call check(maxval(abs(r%x-[0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]))<2.0e-5_dp, &
        'six-binary solution')
    print '(a)', 'PASS test_mixed_integer'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; print '(a,1x,a)','FAIL',msg; error stop 1; end if
    end subroutine
end program
