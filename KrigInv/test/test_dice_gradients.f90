! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! Vendored in KrigInv-fortran under the GPL-3 option.
program test_dice_gradients
  use dk_kinds, only : dp
  use dk_covariance
  use dk_model
  implicit none
  integer, parameter :: n=10, d=2
  real(dp) :: x(n,d), y(n), q(3), g(3), gf(3), p(2), gl(2), glf(2), h, fp, fm
  real(dp), allocatable :: f(:,:), noise(:)
  type(km_model) :: m
  type(covariance_model) :: cov
  integer :: i, k

  do i=1,n
    x(i,1)=mod(real(7*i,dp)*0.113_dp,1.0_dp)
    x(i,2)=mod(real(11*i,dp)*0.071_dp,1.0_dp)
    y(i)=cos(1.3_dp*x(i,1))+0.2_dp*x(i,2)+0.05_dp*sin(real(i,dp))
  end do
  call trend_linear(x,f)
  cov%kind=covariance_kind('gauss'); allocate(cov%range(d)); cov%range=[0.3_dp,0.5_dp]; cov%sd2=0.8_dp
  allocate(noise(n)); noise=0.03_dp
  call km_fit_cov(m,x,y,f,cov,estimate_cov=.false.,estimate_var=.false.,estimate_trend=.true.,noise_var=noise)

  m%estimate_cov=.true.; m%estimate_var=.true.; m%fit_case=3
  q=[0.3_dp,0.5_dp,0.8_dp]
  call loglik_grad(m,q,g)
  do k=1,3
    h=1.0e-6_dp
    q(k)=q(k)+h; fp=loglik_fun(m,q)
    q(k)=q(k)-2.0_dp*h; fm=loglik_fun(m,q)
    q(k)=q(k)+h; gf(k)=(fp-fm)/(2.0_dp*h)
  end do
  if(maxval(abs(g-gf))>5.0e-5_dp) error stop 'noisy likelihood gradient mismatch'

  if(allocated(m%noise_var)) deallocate(m%noise_var)
  m%noise_flag=.false.; m%fit_case=1; p=[0.3_dp,0.5_dp]
  call leave_one_out_grad(m,p,gl)
  do k=1,2
    h=1.0e-6_dp
    p(k)=p(k)+h; fp=leave_one_out_fun(m,p)
    p(k)=p(k)-2.0_dp*h; fm=leave_one_out_fun(m,p)
    p(k)=p(k)+h; glf(k)=(fp-fm)/(2.0_dp*h)
  end do
  if(maxval(abs(gl-glf))>5.0e-5_dp) error stop 'LOO gradient mismatch'

  print *, 'test_gradients: PASS'
end program test_dice_gradients
