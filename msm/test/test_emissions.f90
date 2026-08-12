program test_emissions
    use msm, only : dp, emission_model, emission_probability, emission_derivative, &
        hmm_cat,hmm_ident,hmm_unif,hmm_norm,hmm_lnorm,hmm_exp,hmm_gamma,hmm_weibull,hmm_pois,hmm_binom, &
        hmm_beta_binom,hmm_tnorm,hmm_metnorm,hmm_meunif,hmm_nbinom,hmm_beta,hmm_t
    implicit none
    type(emission_model) :: m,mp,mm
    real(dp),allocatable :: d(:)
    real(dp) :: f,fd,eps,x
    integer :: kinds(17),i,j
    kinds=[hmm_cat,hmm_ident,hmm_unif,hmm_norm,hmm_lnorm,hmm_exp,hmm_gamma,hmm_weibull,hmm_pois,hmm_binom, &
        hmm_beta_binom,hmm_tnorm,hmm_metnorm,hmm_meunif,hmm_nbinom,hmm_beta,hmm_t]
    m%kind=hmm_norm; m%pars=[1.2_dp,0.7_dp]; x=0.8_dp; d=emission_derivative(m,x); eps=1e-6_dp
    do j=1,2
        mp=m; mm=m; mp%pars(j)=mp%pars(j)+eps; mm%pars(j)=mm%pars(j)-eps
        fd=(emission_probability(mp,x)-emission_probability(mm,x))/(2.0_dp*eps)
        call check(abs(d(j)-fd)<2e-8_dp,"normal derivative")
    end do
    do i=1,17
        call setup(kinds(i),m,x)
        f=emission_probability(m,x)
        call check(f>=0.0_dp .and. f<huge(1.0_dp),"finite emission density")
    end do
    m%kind=hmm_cat; m%pars=[0.2_dp,0.3_dp,0.5_dp]
    call check(abs(emission_probability(m,3.0_dp)-0.5_dp)<1e-15_dp,"categorical")
    print '(a)', 'test_emissions: PASS'
contains
    subroutine setup(kind,m,x)
        integer,intent(in)::kind; type(emission_model),intent(out)::m; real(dp),intent(out)::x
        m%kind=kind; x=1.0_dp
        select case(kind)
        case(hmm_cat); m%pars=[0.3_dp,0.7_dp]
        case(hmm_ident); m%pars=[1.0_dp]
        case(hmm_unif); m%pars=[0.0_dp,2.0_dp]
        case(hmm_norm); m%pars=[0.0_dp,1.2_dp]
        case(hmm_lnorm); m%pars=[0.0_dp,0.8_dp]
        case(hmm_exp); m%pars=[1.4_dp]
        case(hmm_gamma); m%pars=[2.0_dp,1.3_dp]
        case(hmm_weibull); m%pars=[1.5_dp,2.0_dp]
        case(hmm_pois); m%pars=[2.5_dp]
        case(hmm_binom); m%pars=[5.0_dp,0.3_dp]
        case(hmm_beta_binom); m%pars=[5.0_dp,0.4_dp,0.2_dp]
        case(hmm_tnorm); m%pars=[0.0_dp,1.0_dp,-2.0_dp,2.0_dp]
        case(hmm_metnorm); m%pars=[0.0_dp,1.0_dp,-2.0_dp,2.0_dp,0.5_dp,0.0_dp]
        case(hmm_meunif); m%pars=[0.0_dp,2.0_dp,0.5_dp,0.0_dp]
        case(hmm_nbinom); m%pars=[3.0_dp,0.6_dp]
        case(hmm_beta); m%pars=[2.0_dp,3.0_dp]; x=0.4_dp
        case(hmm_t); m%pars=[0.0_dp,1.5_dp,5.0_dp]
        end select
    end subroutine setup
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; write(*,'(a)') 'FAIL: '//msg; error stop 1; end if
    end subroutine check
end program test_emissions
