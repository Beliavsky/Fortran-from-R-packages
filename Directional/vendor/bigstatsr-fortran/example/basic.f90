program basic
    use bigstatsr
    implicit none
    type(fbm_real), target :: x
    type(colstats_result) :: st
    type(big_svd_result) :: fit
    real(dp) :: a(5,3)
    a=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp, &
               2.0_dp,1.0_dp,0.0_dp,-1.0_dp,-2.0_dp, &
               1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp],[5,3])
    x=create_fbm('basic_bigstatsr.bk',5,3)
    call fbm_from_array(x,a)
    st=big_colstats(x)
    fit=big_svd(x,2,center=.true.)
    print '(a,3f10.4)', 'column means: ',st%sum/5.0_dp
    print '(a,2f10.4)', 'singular values: ',fit%d
    call execute_command_line('rm -f basic_bigstatsr.bk')
end program basic
