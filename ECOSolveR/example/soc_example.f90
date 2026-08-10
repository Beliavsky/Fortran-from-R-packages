program soc_example
    use ecos
    implicit none
    type(ecos_dims)::dims
    type(ecos_result)::res
    real(dp)::c(3),g(3,3),h(3),a(2,3),b(2)
    allocate(dims%q(1)); dims%q=[3]
    c=[1.0_dp,0.0_dp,0.0_dp]
    g=0.0_dp; g(1,1)=-1.0_dp; g(2,2)=-1.0_dp; g(3,3)=-1.0_dp
    h=0.0_dp; a=0.0_dp; a(1,2)=1.0_dp; a(2,3)=1.0_dp; b=[3.0_dp,4.0_dp]
    call ecos_csolve(c,g,h,dims,res,a,b)
    print '(a,i0)', 'status: ',res%exitflag
    print '(a,3f14.8)', 'x: ',res%x
    print '(a,f14.8)', 'objective: ',res%pcost
end program
