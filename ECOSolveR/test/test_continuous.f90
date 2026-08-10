program test_continuous
    use ecos
    implicit none
    type(ecos_dims) :: dims
    type(ecos_result) :: res
    real(dp) :: c2(2), g2(2,2), h2(2), a1(1,2), b1(1)
    real(dp) :: c3(3), g3(3,3), h3(3), a2(2,3), b2(2)

    c2 = [-1.0_dp,-1.0_dp]
    g2 = 0.0_dp
    g2(1,1)=-1.0_dp; g2(2,2)=-1.0_dp
    h2=0.0_dp; a1=1.0_dp; b1=1.0_dp
    dims%l=2
    call ecos_csolve(c2,g2,h2,dims,res,a1,b1)
    call check(res%exitflag==ECOS_OPTIMAL,'LP status')
    call check(abs(sum(res%x)-1.0_dp)<1.0e-7_dp,'LP equality')
    call check(abs(res%pcost+1.0_dp)<1.0e-7_dp,'LP objective')

    dims%l=0; dims%e=0
    allocate(dims%q(1)); dims%q=[3]
    c3=[1.0_dp,0.0_dp,0.0_dp]
    g3=0.0_dp; g3(1,1)=-1.0_dp; g3(2,2)=-1.0_dp; g3(3,3)=-1.0_dp
    h3=0.0_dp; a2=0.0_dp; a2(1,2)=1.0_dp; a2(2,3)=1.0_dp; b2=[3.0_dp,4.0_dp]
    call ecos_csolve(c3,g3,h3,dims,res,a2,b2)
    call check(res%exitflag==ECOS_OPTIMAL,'SOC status')
    call check(maxval(abs(res%x-[5.0_dp,3.0_dp,4.0_dp]))<2.0e-7_dp,'SOC solution')

    deallocate(dims%q); dims%e=1
    c3=[0.0_dp,1.0_dp,0.0_dp]
    a2=0.0_dp; a2(1,1)=1.0_dp; a2(2,3)=1.0_dp; b2=[1.0_dp,1.0_dp]
    call ecos_csolve(c3,g3,h3,dims,res,a2,b2, &
        control=ecos_control(maxit=150,feastol=1.0e-9_dp,reltol=1.0e-9_dp,abstol=1.0e-9_dp))
    call check(res%exitflag==ECOS_OPTIMAL,'exponential status')
    call check(abs(res%x(2)-exp(1.0_dp))<1.0e-7_dp,'exponential solution')

    print '(a)', 'PASS test_continuous'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(*),intent(in)::msg
        if(.not.ok) then; print '(a,1x,a)','FAIL',msg; error stop 1; end if
    end subroutine
end program
