program test_hmm_extended
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use msm, only : dp, emission_model, hmm_cat, hmm_norm, make_generator, transition_matrix, &
        hmm_result, hmm_forward_backward, obs_panel, obs_death
    implicit none
    type(emission_model) :: mod(2,2),dm(1,3)
    type(hmm_result) :: r
    real(dp) :: off(2,2),q3(2,2,1),tm(3),ob(2,3),ip(2),brute,pr
    real(dp) :: offd(3,3),qd(3,3,1),td(2),od(1,2),ipd(3),expected
    real(dp),allocatable :: p(:,:)
    logical :: allowed(2,3)
    integer :: ot(2),ds(2),a,b,c

    off=0.0_dp; off(1,2)=0.25_dp; off(2,1)=0.10_dp; q3(:,:,1)=make_generator(off)
    mod(1,1)%kind=hmm_cat; mod(1,1)%pars=[0.85_dp,0.15_dp]
    mod(1,2)%kind=hmm_cat; mod(1,2)%pars=[0.20_dp,0.80_dp]
    mod(2,1)%kind=hmm_norm; mod(2,1)%pars=[0.0_dp,1.0_dp]
    mod(2,2)%kind=hmm_norm; mod(2,2)%pars=[2.0_dp,1.0_dp]
    tm=[0.0_dp,1.0_dp,2.0_dp]; ob(:,1)=[1.0_dp,0.2_dp]; ob(:,2)=[2.0_dp,ieee_value(0.0_dp,ieee_quiet_nan)]; ob(:,3)=[2.0_dp,1.8_dp]
    ip=[0.6_dp,0.4_dp]; allowed=.true.; allowed(2,2)=.false.
    call hmm_forward_backward(q3,tm,ob,mod,ip,r,allowed_state=allowed)
    brute=0.0_dp
    do a=1,2; do b=1,2; do c=1,2
        if(.not.allowed(b,2)) cycle
        p=transition_matrix(q3(:,:,1),1.0_dp)
        pr=ip(a)*emit(a,1)*p(a,b)*emit(b,2)*p(b,c)*emit(c,3)
        brute=brute+pr
    end do; end do; end do
    call check(abs(r%minus2loglik+2.0_dp*log(brute))<5e-13_dp,"multivariate/missing/allowed HMM")

    offd=0.0_dp; offd(1,2)=0.2_dp; offd(1,3)=0.1_dp; offd(2,3)=0.4_dp; qd(:,:,1)=make_generator(offd)
    dm(1,1)%kind=hmm_cat; dm(1,1)%pars=[1.0_dp]
    dm(1,2)%kind=hmm_cat; dm(1,2)%pars=[1.0_dp]
    dm(1,3)%kind=hmm_cat; dm(1,3)%pars=[1.0_dp]
    td=[0.0_dp,2.0_dp]; od=1.0_dp; ipd=[1.0_dp,0.0_dp,0.0_dp]; ot=[obs_panel,obs_death]; ds=[0,3]
    call hmm_forward_backward(qd,td,od,dm,ipd,r,obstype=ot,death_state=ds)
    p=transition_matrix(qd(:,:,1),2.0_dp)
    expected=p(1,1)*qd(1,3,1)+p(1,2)*qd(2,3,1)
    call check(abs(r%minus2loglik+2.0_dp*log(expected))<5e-13_dp,"exact-death HMM")
    print '(a)', 'test_hmm_extended: PASS'
contains
    function emit(st,k) result(v)
        integer,intent(in)::st,k
        real(dp)::v,z,mean,sd
        integer::cat
        cat=nint(ob(1,k)); v=mod(1,st)%pars(cat)
        if(.not.ieee_is_nan(ob(2,k))) then
            mean=mod(2,st)%pars(1); sd=mod(2,st)%pars(2); z=(ob(2,k)-mean)/sd
            v=v*exp(-0.5_dp*z*z)/(sd*sqrt(2.0_dp*acos(-1.0_dp)))
        end if
    end function emit
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; write(*,'(a)') 'FAIL: '//msg; error stop 1; end if
    end subroutine check
end program test_hmm_extended
