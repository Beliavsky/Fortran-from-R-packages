program hidden_markov
    use msm, only : dp, make_generator, emission_model, hmm_cat, hmm_result, hmm_forward_backward
    implicit none
    real(dp) :: off(2,2), q3(2,2,1), times(5), obs(1,5), initp(2)
    type(emission_model) :: models(1,2)
    type(hmm_result) :: fit
    integer :: k

    off = 0.0_dp
    off(1,2) = 0.25_dp
    off(2,1) = 0.10_dp
    q3(:,:,1) = make_generator(off)

    models(1,1)%kind = hmm_cat
    models(1,1)%pars = [0.90_dp, 0.10_dp]
    models(1,2)%kind = hmm_cat
    models(1,2)%pars = [0.15_dp, 0.85_dp]

    times = [0.0_dp, 1.0_dp, 2.5_dp, 4.0_dp, 6.0_dp]
    obs(1,:) = [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp]
    initp = [0.8_dp, 0.2_dp]

    call hmm_forward_backward(q3, times, obs, models, initp, fit)
    print '(a,f12.6)', '-2 log likelihood = ', fit%minus2loglik
    print '(a,*(i0,1x))', 'Viterbi states: ', fit%viterbi
    print '(a)', 'Smoothed state-2 probabilities:'
    do k = 1, size(times)
        print '(f6.2,1x,f10.6)', times(k), fit%smoothed(2,k)
    end do
end program hidden_markov
