program sample_polytope
    use limsolve
    implicit none
    type(sample_result) :: s
    real(dp) :: a0(0,2),b0(0),e(1,2),f(1),g(2,2),h(2)
    e(1,:)=[1.0_dp,1.0_dp]; f=[1.0_dp]
    g=0.0_dp; g(1,1)=1.0_dp; g(2,2)=1.0_dp; h=0.0_dp
    call xsample(a0,b0,e,f,g,h,s,iter=1000,outputlength=100,type='rda',seed=1234)
    print '(a,f10.5)', 'accepted ratio = ',s%accepted_ratio
    print '(a,2f12.6)', 'sample mean = ',sum(s%x(:,1))/size(s%x,1), &
        sum(s%x(:,2))/size(s%x,1)
end program sample_polytope
