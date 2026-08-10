program mixed_integer_example
    use ecos
    implicit none
    type(ecos_dims)::dims
    type(ecos_result)::res
    real(dp)::c(2),g(2,2),h(2)
    integer::bool_vars(1)
    dims%l=2
    c=[-1.1_dp,-1.0_dp]
    g=reshape([2.0_dp,3.0_dp,1.0_dp,4.0_dp],[2,2]); h=[4.0_dp,12.0_dp]
    bool_vars=[1]
    call ecos_csolve(c,g,h,dims,res,bool_vars=bool_vars)
    print '(a,i0)', 'status: ',res%exitflag
    print '(a,2f14.8)', 'x: ',res%x
    print '(a,f14.8)', 'objective: ',res%pcost
    print '(a,i0)', 'B&B nodes: ',res%mi_iter
end program
