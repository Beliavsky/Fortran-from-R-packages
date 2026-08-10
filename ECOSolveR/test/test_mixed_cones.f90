program test_mixed_cones
    use ecos
    implicit none
    type(ecos_dims)::d
    type(ecos_result)::r
    real(dp)::c(4),g(4,4),h(4),a(2,4),b(2)
    ! x1 >= 0 plus (t,x2,x3) in Q3, with x2=6,x3=8 and x1=t.
    d%l=1; allocate(d%q(1)); d%q=[3]
    c=[1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    g=0.0_dp; g(1,1)=-1.0_dp; g(2,1)=-1.0_dp; g(3,2)=-1.0_dp; g(4,3)=-1.0_dp
    h=0.0_dp
    a=0.0_dp; a(1,2)=1.0_dp; a(2,3)=1.0_dp; b=[6.0_dp,8.0_dp]
    call ecos_csolve(c,g,h,d,r,a,b)
    call check(r%exitflag==ECOS_OPTIMAL,'mixed cone status')
    call check(abs(r%x(1)-10.0_dp)<1.0e-6_dp,'mixed cone norm')
    print '(a)', 'PASS test_mixed_cones'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; print '(a,1x,a)','FAIL',msg; error stop 1; end if
    end subroutine
end program
