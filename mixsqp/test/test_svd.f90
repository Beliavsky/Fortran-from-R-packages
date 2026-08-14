program test_svd
  use mixsqp
  implicit none
  integer,parameter::n=60,m=8,rk=3
  real(dp)::A(n,rk),B(m,rk),L(n,m)
  type(mixsqp_result)::full,low
  type(mixsqp_control)::c1,c2
  integer::i,j,k
  do k=1,rk
    do i=1,n
      A(i,k)=0.2_dp+sin(0.07_dp*real(i*k,dp))**2
    end do
    do j=1,m
      B(j,k)=0.1_dp+cos(0.11_dp*real(j*(k+1),dp))**2
    end do
  end do
  L=matmul(A,transpose(B))
  c1=mixsqp_default_control();c1%tol_svd=0._dp;c1%verbose=.false.
  c2=c1;c2%tol_svd=1e-10_dp
  call fit_mixsqp(L,full,control=c1)
  call fit_mixsqp(L,low,control=c2)
  call check(full%status==0 .and. low%status==0,'SVD test did not converge')
  call check(low%used_svd,'low-rank path not used')
  call check(low%svd_rank==rk,'wrong detected rank')
  call check(abs(full%value-low%value)<1e-9_dp,'SVD objective mismatch')
  call check(maxval(abs(full%x-low%x))<2e-6_dp,'SVD solution mismatch')
  print *, 'test_svd: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) error stop msg
  end subroutine
end program
