
program test_joker
  use joker
  implicit none
  integer :: fails, i
  real(dp), allocatable :: xs(:), dm(:,:), mg(:,:), draw(:)
  integer, allocatable :: xi(:)
  type(estimate2) :: e2
  type(estimate1) :: e1
  type(estimate_vec) :: ev
  type(estimate_vec_scale) :: emg
  real(dp) :: alpha(3), point(3), shape(3)
  fails=0

  call check_close(pnorm_j(0.0_dp,0.0_dp,1.0_dp),0.5_dp,1e-12_dp,"normal cdf")
  call check_close(qnorm_j(0.975_dp,0.0_dp,1.0_dp),1.95996398454005_dp,2e-7_dp,"normal quantile")
  call check_close(pbeta(0.5_dp,2.0_dp,3.0_dp),0.6875_dp,1e-10_dp,"beta cdf")
  call check_close(pgamma_j(qgamma_j(.37_dp,2.5_dp,1.7_dp),2.5_dp,1.7_dp),.37_dp,2e-10_dp,"gamma inversion")
  call check_close(pt_j(qt_j(.8_dp,7.0_dp),7.0_dp),.8_dp,2e-10_dp,"student inversion")
  call check_close(pfisher(qfisher(.6_dp,5.0_dp,10.0_dp),5.0_dp,10.0_dp),.6_dp,2e-10_dp,"F inversion")
  call check_close(plaplace(qlaplace(.23_dp,2.0_dp,.7_dp),2.0_dp,.7_dp),.23_dp,1e-12_dp,"laplace inversion")

  allocate(xs(5)); xs=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
  e2=enorm(xs)
  call check_close(e2%p1,3.0_dp,1e-12_dp,"normal MLE mean")
  call check_close(e2%p2,sqrt(2.0_dp),1e-12_dp,"normal MLE sd")
  e1=eexp(xs)
  call check_close(e1%p1,1.0_dp/3.0_dp,1e-12_dp,"exp MLE")

  deallocate(xs); allocate(xs(6)); xs=[.1_dp,.2_dp,.3_dp,.4_dp,.6_dp,.7_dp]
  e2=ebeta_mle(xs)
  if(e2%p1<=0.or.e2%p2<=0)call fail("beta mle positive")
  e2=egamma_mle([.5_dp,1.0_dp,1.5_dp,2.0_dp,3.0_dp,4.0_dp])
  if(e2%p1<=0.or.e2%p2<=0)call fail("gamma mle positive")
  e2=eweibull_mle([.5_dp,.8_dp,1.1_dp,1.7_dp,2.4_dp,3.0_dp])
  if(e2%p1<=0.or.e2%p2<=0)call fail("weibull mle positive")

  alpha=[2.0_dp,3.0_dp,4.0_dp]
  point=[.2_dp,.3_dp,.5_dp]
  if(ddir(point,alpha)<=0)call fail("dirichlet density")
  allocate(dm(5000,3),draw(3))
  do i=1,size(dm,1)
    call rdir(alpha,draw); dm(i,:)=draw
  end do
  call check_close(sum(dm(:,1))/size(dm,1),2.0_dp/9.0_dp,.02_dp,"dirichlet rng mean")
  ev=edir_me(dm)
  if(any(ev%value<=0))call fail("dirichlet ME positive")
  ev=edir_mle(dm)
  if(any(ev%value<=0))call fail("dirichlet MLE positive")

  shape=[1.5_dp,2.0_dp,3.0_dp]
  allocate(mg(4000,3))
  do i=1,size(mg,1)
    call rmultigam(shape,2.0_dp,draw); mg(i,:)=draw
  end do
  call check_close(sum(mg(:,3))/size(mg,1),sum(shape)*2.0_dp,.35_dp,"multigamma rng mean")
  emg=emultigam_me(mg)
  if(any(emg%value<=0).or.emg%scale<=0)call fail("multigamma ME positive")
  emg=emultigam_mle(mg)
  if(any(emg%value<=0).or.emg%scale<=0)call fail("multigamma MLE positive")

  allocate(xi(6)); xi=[0,1,1,0,1,1]
  e1=ebern(xi); call check_close(e1%p1,2.0_dp/3.0_dp,1e-12_dp,"bernoulli estimator")

  if(fails==0)then
    print '(a)',"test_joker: PASS"
  else
    print '(a,i0)',"test_joker: FAIL ",fails
    error stop 1
  end if
contains
  subroutine check_close(a,b,tol,name)
    real(dp),intent(in)::a,b,tol
    character(*),intent(in)::name
    if(abs(a-b)>tol)then
      print '(a,2es16.7)',trim(name)//" mismatch: ",a,b
      fails=fails+1
    end if
  end subroutine
  subroutine fail(name)
    character(*),intent(in)::name
    print '(a)',trim(name)
    fails=fails+1
  end subroutine
end program
