program test_core
  use anmc
  implicit none
  real(dp), allocatable :: x(:,:), xt(:,:), e(:,:), sigma(:,:), mu(:), pn(:)
  real(dp) :: m1, m2, v1, v2, c12, exact_half_mean
  integer, allocatable :: idx(:)
  type(active_dims_result) :: aq
  type(probability_estimate) :: pmax, pmin
  type(conservative_result) :: ce
  type(anmc_problem) :: problem
  type(mc_params) :: pars
  type(mc_result) :: mr, ar
  type(simulation_control) :: sctl
  logical :: ok
  integer :: i

  call seed_fortran_rng(12345)
  mu = [1.0_dp,-2.0_dp]
  sigma = reshape([1.0_dp,0.3_dp,0.3_dp,2.0_dp],[2,2])
  x = mvrnorm_arma(30000,mu,sigma,0,ok)
  call assert_true(ok,'mvrnorm_arma failed')
  m1=sum(x(1,:))/real(size(x,2),dp); m2=sum(x(2,:))/real(size(x,2),dp)
  v1=sum((x(1,:)-m1)**2)/real(size(x,2)-1,dp)
  v2=sum((x(2,:)-m2)**2)/real(size(x,2)-1,dp)
  c12=sum((x(1,:)-m1)*(x(2,:)-m2))/real(size(x,2)-1,dp)
  call assert_close(m1,1.0_dp,0.03_dp,'normal mean 1')
  call assert_close(m2,-2.0_dp,0.04_dp,'normal mean 2')
  call assert_close(v1,1.0_dp,0.05_dp,'normal variance 1')
  call assert_close(v2,2.0_dp,0.08_dp,'normal variance 2')
  call assert_close(c12,0.3_dp,0.05_dp,'normal covariance')
  call seed_fortran_rng(12345)
  x=mvrnorm_arma(30000,mu,reshape([1.0_dp,0.0_dp,0.3_dp,sqrt(1.91_dp)],[2,2]),1,ok)
  call assert_true(ok,'mvrnorm_arma Cholesky-input mode failed')
  m1=sum(x(1,:))/real(size(x,2),dp); m2=sum(x(2,:))/real(size(x,2),dp)
  c12=sum((x(1,:)-m1)*(x(2,:)-m2))/real(size(x,2)-1,dp)
  call assert_close(c12,0.3_dp,0.05_dp,'Cholesky-input covariance')

  call seed_fortran_rng(42)
  xt=trmvrnorm_rej_cpp(6000,[0.0_dp],reshape([1.0_dp],[1,1]),[0.0_dp],[huge(1.0_dp)],0,ok=ok)
  call assert_true(ok,'truncated normal failed')
  exact_half_mean=sqrt(2.0_dp/acos(-1.0_dp))
  call assert_close(sum(xt(1,:))/real(size(xt,2),dp),exact_half_mean,0.04_dp,'truncated normal mean')
  call assert_true(minval(xt)>=0.0_dp,'truncated sample violates lower bound')

  allocate(e(10,1)); e(:,1)=[(real(i,dp),i=1,10)]
  mu=[(0.0_dp,i=1,10)]; sigma=0.5_dp
  deallocate(sigma); allocate(sigma(10,10)); sigma=0.0_dp
  do i=1,10; sigma(i,i)=1.0_dp; end do
  idx=select_active_dims(3,e,0.0_dp,mu,sigma,method=0)
  call assert_true(all(idx==[1,5,10]),'method-0 active dimensions')
  aq=select_q_dims(e,0.0_dp,mu,sigma,method=0,limits=[4,4],reduced_return=.false.)
  call assert_true(size(aq%ind_q)==4,'select_q_dims fixed range')
  call assert_true(all(aq%ind_q==[1,4,7,10]),'select_q_dims indices')

  deallocate(mu,sigma,e)
  allocate(mu(4),sigma(4,4),e(4,1)); mu=0.0_dp; e(:,1)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp]
  sigma=0.5_dp
  do i=1,4; sigma(i,i)=1.0_dp; end do
  pmax=proba_max(0.01_dp,0.0_dp,mu,sigma,e,q=4,method=0,prob_control=genz_bretz(maxpts=200000,abseps=1e-5_dp))
  pmin=proba_min(0.01_dp,0.0_dp,mu,sigma,e,q=4,method=0,prob_control=genz_bretz(maxpts=200000,abseps=1e-5_dp))
  call assert_true(pmax%ok .and. pmin%ok,'orthant probability call failed')
  call assert_close(pmax%probability,0.8_dp,0.01_dp,'equicorrelated max probability')
  call assert_close(pmin%probability,0.8_dp,0.01_dp,'equicorrelated min probability')

  problem%mu_eq=[0.0_dp]
  problem%sigma_eq=reshape([1.0_dp],[1,1])
  problem%threshold=0.0_dp
  problem%mu_emq=[0.0_dp]
  problem%ww_cond_q=reshape([0.0_dp],[1,1])
  problem%sigma_cond_q_chol=reshape([1.0_dp],[1,1])
  pars%cx=1.0e-5_dp; pars%beta=1.0e-5_dp; pars%alpha_cost=0.0_dp; pars%eval_g=0.0_dp
  sctl%max_outer=5000; sctl%max_inner=100; sctl%enforce_budget=.false.
  call seed_fortran_rng(777)
  mr=mc_gauss(0.02_dp,problem,params=pars,sim_control=sctl)
  call assert_true(mr%ok,'MC_Gauss failed')
  call assert_close(mr%estim,0.5_dp,0.06_dp,'MC conditional probability')
  call seed_fortran_rng(777)
  ar=anmc_gauss(0.02_dp,problem,sim_control=sctl,fixed_n=2500,fixed_m=20)
  call assert_true(ar%ok,'ANMC_Gauss failed')
  call assert_close(ar%estim,0.5_dp,0.06_dp,'ANMC conditional probability')
  call assert_true(ar%var_est>=0.0_dp,'ANMC variance negative')

  pn=[0.99_dp,0.95_dp,0.8_dp,0.3_dp]
  ce=conservative_estimate(0.8_dp,mu,sigma,e,0.0_dp,pn=pn,excursion_type='>',prob_control=genz_bretz(maxpts=100000))
  call assert_true(ce%ok,'conservative estimate failed')
  call assert_true(size(ce%set)==4,'conservative set shape')
  call assert_true(ce%level>=0.0_dp .and. ce%level<=1.0_dp,'conservative level bounds')
  pn=[0.99_dp,0.98_dp,0.97_dp,0.96_dp]
  ce=conservative_estimate(0.90_dp,mu,sigma,e,0.0_dp,pn=pn,excursion_type='>', &
                           prob_control=genz_bretz(maxpts=100000))
  call assert_true(ce%ok,'all-high-coverage conservative estimate failed')
  call assert_true(size(ce%set)==4,'all-high-coverage set shape')

  print *, 'test_core: PASS'
contains
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    if(abs(x-y)>tol) then
      write(*,*) 'FAIL: ',trim(msg),' got=',x,' expected=',y,' tol=',tol
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(cond,msg)
    logical,intent(in)::cond
    character(len=*),intent(in)::msg
    if(.not.cond) then
      write(*,*) 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine assert_true
end program test_core
