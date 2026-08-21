program bootstrap_mean
    use boot_kinds, only : dp
    use boot_core
    use boot_statistics, only : mean_dp, variance_dp
    use boot_ci, only : percentile_ci
    implicit none
    real(dp)::data(8,1),lo(1),hi(1)
    type(bootstrap_result)::res
    data(:,1)=[2.0_dp,3.0_dp,5.0_dp,7.0_dp,11.0_dp,13.0_dp,17.0_dp,19.0_dp]
    call bootstrap_weighted(data,stat,999,'ordinary',res)
    call percentile_ci(res%t,[0.95_dp],lo,hi)
    print '(a,f10.5)', 'observed mean = ',res%t0
    print '(a,f10.5)', 'bootstrap se  = ',sqrt(variance_dp(res%t))
    print '(a,2f10.5)', '95% percentile CI = ',lo(1),hi(1)
contains
    function stat(x,w) result(v)
        real(dp),intent(in)::x(:,:),w(:)
        real(dp)::v
        v=sum(x(:,1)*w)/sum(w)
    end function stat
end program bootstrap_mean
