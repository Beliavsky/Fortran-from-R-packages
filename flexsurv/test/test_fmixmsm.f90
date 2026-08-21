program test_fmixmsm
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_spec, initialize_spec, initial_theta
  use flexsurv_distributions, only : dist_exponential
  use flexsurv_mixture_full, only : flexsurvmix_full_result
  use flexsurv_fmixmsm, only : fmix_transition_model, fmixmsm_model, build_fmixmsm, &
    fmix_path_probabilities, fmix_path_mean_times, fmix_final_probabilities, &
    fmix_final_mean_times, fmix_path_time_quantiles, &
    fmix_path_probabilities_ci, fmix_final_mean_times_ci
  implicit none
  type(fmix_transition_model) :: node(2)
  type(fmixmsm_model) :: msm
  real(dp) :: pp(2),mp(2),fp(3),fm(3),qq(2,3),lo2(2),hi2(2),lo3(3),hi3(3)

  node(1)%from=1;allocate(node(1)%to(2));node(1)%to=[2,3]
  call make_mix(node(1)%mix,[0.4_dp,0.6_dp],[0.5_dp,0.25_dp])
  node(2)%from=2;allocate(node(2)%to(1));node(2)%to=[3]
  call make_mix(node(2)%mix,[1.0_dp],[1.0_dp])
  msm=build_fmixmsm(node,3,1)
  if(msm%pathways%has_cycle) error stop 'unexpected cycle'
  if(msm%pathways%npath/=2) error stop 'path count'
  call fmix_path_probabilities(msm,1,pp)
  call fmix_path_mean_times(msm,1,mp)
  if(abs(sum(pp)-1.0_dp)>1.0e-12_dp) error stop 'path prob sum'
  if(maxval(abs(sort2(pp)-[0.4_dp,0.6_dp]))>1.0e-10_dp) error stop 'path probs'
  if(maxval(abs(sort2(mp)-[3.0_dp,4.0_dp]))>1.0e-10_dp) error stop 'path means'
  call fmix_final_probabilities(msm,1,fp)
  call fmix_final_mean_times(msm,1,fm)
  if(abs(fp(3)-1.0_dp)>1.0e-12_dp) error stop 'final probability'
  if(abs(fm(3)-3.6_dp)>1.0e-10_dp) error stop 'final mean'
  call fmix_path_time_quantiles(msm,1,1000,[0.25_dp,0.5_dp,0.75_dp],qq,123)
  if(any(qq<=0.0_dp)) error stop 'path quantiles'
  call fmix_path_probabilities_ci(msm,1,20,0.95_dp,pp,lo2,hi2,321)
  if(any(lo2>pp).or.any(hi2<pp)) error stop 'path bootstrap ci'
  call fmix_final_mean_times_ci(msm,1,20,0.95_dp,fm,lo3,hi3,654)
  if(lo3(3)>fm(3).or.hi3(3)<fm(3)) error stop 'mean bootstrap ci'
  print *, 'test_fmixmsm: PASS'
contains
  subroutine make_mix(mix,prob,rate)
    type(flexsurvmix_full_result),intent(out)::mix
    real(dp),intent(in)::prob(:),rate(:)
    type(flexsurv_spec)::sp
    integer::k,j
    k=size(prob);mix%k=k;mix%nobs=1;mix%nprob_cov=0
    allocate(mix%specs(k),mix%components(k),mix%alpha(max(0,k-1)),mix%prob_beta(0,max(0,k-1)),mix%prob_x(1,0))
    do j=1,k
      call initialize_spec(sp,dist_exponential,1,[rate(j)]);mix%specs(j)=sp
      mix%components(j)%theta=initial_theta(sp)
    end do
    if(k>1)then
      do j=2,k;mix%alpha(j-1)=log(prob(j)/prob(1));end do
    end if
    mix%npar=max(0,k-1)+k
    allocate(mix%covariance(mix%npar,mix%npar));mix%covariance=0.0_dp
    do j=1,mix%npar;mix%covariance(j,j)=0.0025_dp;end do
  end subroutine make_mix
  pure function sort2(x) result(y)
    real(dp),intent(in)::x(2);real(dp)::y(2)
    if(x(1)<=x(2))then;y=x;else;y=[x(2),x(1)];end if
  end function sort2
end program test_fmixmsm
