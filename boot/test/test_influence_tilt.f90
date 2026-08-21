program test_influence_tilt
    use boot_kinds, only : dp
    use boot_influence
    use boot_tilt
    implicit none
    real(dp)::data(5,1),l(5),lj(5),lp(5),p(5),lambda,target
    integer::strata(5),info
    data(:,1)=[1.0_dp,2.0_dp,4.0_dp,7.0_dp,11.0_dp]
    strata=1
    call infinitesimal_jackknife(data,stat,strata,l)
    call delete1_jackknife(data,stat,strata,lj)
    call positive_jackknife(data,stat,strata,lp)
    if(maxval(abs(l-(data(:,1)-5.0_dp)))>2.0e-3_dp)error stop 1
    if(maxval(abs(lj-(data(:,1)-5.0_dp)))>1.0e-12_dp)error stop 2
    if(maxval(abs(lp-(data(:,1)-5.0_dp)))>1.0e-12_dp)error stop 3
    target=0.4_dp
    call exponential_tilt(l,target,0.0_dp,strata,p,lambda,info)
    if(info/=0)error stop 4
    if(abs(sum(l*p)-target)>1.0e-8_dp)error stop 5
    if(abs(sum(p)-1.0_dp)>1.0e-12_dp)error stop 6
    print '(a)', 'test_influence_tilt: PASS'
contains
    function stat(x,w) result(v)
        real(dp),intent(in)::x(:,:),w(:)
        real(dp)::v
        v=sum(x(:,1)*w)/sum(w)
    end function stat
end program test_influence_tilt
