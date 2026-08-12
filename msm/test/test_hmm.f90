program test_hmm
    use msm, only : dp, emission_model, hmm_cat, make_generator, transition_matrix, hmm_result, &
        hmm_forward_backward, hmm_minus2loglik, hmm_q_gradient
    implicit none
    type(emission_model) :: mod(1,2)
    type(hmm_result) :: res
    real(dp) :: off(2,2),q(2,2),q3(2,2,1),times(3),obs(1,3),initp(2),brute,pr,p1,eps,fp,fm
    real(dp), allocatable :: p(:,:),dq(:,:,:),g(:)
    integer :: k,st1,st2,st3,bestpath(3)
    real(dp) :: best
    off=0.0_dp; off(1,2)=0.3_dp; off(2,1)=0.15_dp; q=make_generator(off); q3(:,:,1)=q
    times=[0.0_dp,0.8_dp,1.7_dp]; obs(1,:)=[1.0_dp,1.0_dp,2.0_dp]; initp=[0.7_dp,0.3_dp]
    mod(1,1)%kind=hmm_cat; mod(1,1)%pars=[0.9_dp,0.1_dp]
    mod(1,2)%kind=hmm_cat; mod(1,2)%pars=[0.2_dp,0.8_dp]
    call hmm_forward_backward(q3,times,obs,mod,initp,res)
    brute=0.0_dp; best=-1.0_dp; bestpath=1
    do st1=1,2; do st2=1,2; do st3=1,2
        p=transition_matrix(q,times(2)-times(1)); p1=mod(1,st1)%pars(nint(obs(1,1)))
        pr=initp(st1)*p1*p(st1,st2)*mod(1,st2)%pars(nint(obs(1,2)))
        p=transition_matrix(q,times(3)-times(2)); pr=pr*p(st2,st3)*mod(1,st3)%pars(nint(obs(1,3)))
        brute=brute+pr
        if(pr>best) then; best=pr; bestpath=[st1,st2,st3]; end if
    end do; end do; end do
    call check(abs(res%minus2loglik+2.0_dp*log(brute))<3e-13_dp,"forward likelihood")
    call check(all(res%viterbi==bestpath),"Viterbi path")
    do k=1,3; call check(abs(sum(res%filtered(:,k))-1.0_dp)<2e-14_dp,"filtered normalization"); &
        call check(abs(sum(res%smoothed(:,k))-1.0_dp)<2e-14_dp,"smoothed normalization"); end do
    allocate(dq(2,2,1)); dq=0.0_dp; dq(1,1,1)=-1.0_dp; dq(1,2,1)=1.0_dp
    g=hmm_q_gradient(q,dq,times,obs,mod,initp)
    eps=1e-6_dp; off=0.0_dp; off(1,2)=0.3_dp+eps; off(2,1)=0.15_dp; fp=hmm_minus2loglik(make_generator(off),times,obs,mod,initp)
    off(1,2)=0.3_dp-eps; fm=hmm_minus2loglik(make_generator(off),times,obs,mod,initp)
    call check(abs(g(1)-(fp-fm)/(2.0_dp*eps))<3e-7_dp,"HMM q gradient")
    print '(a)', 'test_hmm: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; write(*,'(a)') 'FAIL: '//msg; error stop 1; end if
    end subroutine check
end program test_hmm
