! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
program demo_sn
  use iso_fortran_env, only : int64
  use sn, only : dp, sn_rng_state, rsn, dsn, psn, qsn, selm_result, selm_fit, sn_ok
  implicit none
  type(sn_rng_state) :: rng
  type(selm_result) :: fit
  real(dp), allocatable :: x(:,:),y(:),e(:)
  real(dp) :: p,q
  integer :: i,n,info

  p=psn(0.5_dp,alpha=2.0_dp)
  q=qsn(p,alpha=2.0_dp,info=info)
  write(*,'(a,f12.8)') 'SN density at zero, alpha=2: ',dsn(0.0_dp,alpha=2.0_dp)
  write(*,'(a,f12.8)') 'CDF/quantile round trip:       ',q

  n=100
  allocate(x(n,2),y(n),e(n))
  x(:,1)=1.0_dp
  do i=1,n
    x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
  end do
  call rng%seed(20260804_int64)
  call rsn(rng,e,omega=0.7_dp,alpha=2.2_dp,info=info)
  if(info/=sn_ok) error stop 'random generation failed'
  y=1.0_dp+2.5_dp*x(:,2)+e
  call selm_fit(x,y,'SN',fit,penalty='NONE',max_iter=3000,tol=2.0e-7_dp, &
                start=[1.0_dp,2.5_dp,log(0.7_dp),2.0_dp])
  if(.not.fit%converged) error stop 'fit failed'
  write(*,'(a,2f12.6)') 'Fitted beta:                   ',fit%beta
  write(*,'(a,f12.6)') 'Fitted omega:                  ',fit%omega
  write(*,'(a,f12.6)') 'Fitted alpha:                  ',fit%alpha
  write(*,'(a,f12.6)') 'Log likelihood:                ',fit%log_likelihood
end program demo_sn
