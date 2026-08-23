program test_parity_v02
  use fastmatrix, only: dp, schur_decomp, svd_decomp, matrix_fun, ols_fit_qr, ols_fit_svd, &
    ols_compat_result, mardia_result, harris_result, mardia_test, harris_test
  implicit none
  real(dp) :: a(3,3),t(3,3),q(3,3),wr(3),wi(3),rec(3,3),ident(3,3)
  real(dp) :: ar(2,2),tr(2,2),qr(2,2),wrr(2),wir(2)
  real(dp) :: x(6,2),y(6),utr(3,3),f(3,3)
  real(dp) :: us(6,2),sv(2),vts(2,2),asvd(6,2)
  real(dp) :: xm(8,2)
  type(ols_compat_result) :: rq,rs
  type(mardia_result) :: mr
  type(harris_result) :: hr
  integer :: info,i

  ident=0.0_dp; do i=1,3; ident(i,i)=1.0_dp; end do
  a=reshape([1.0_dp,4.0_dp,7.0_dp, 2.0_dp,5.0_dp,8.0_dp, 3.0_dp,6.0_dp,10.0_dp],[3,3])
  call schur_decomp(a,t,wr,wi,q,info)
  call check(info==0,'schur info')
  rec=matmul(q,matmul(t,transpose(q)))
  call check(maxval(abs(rec-a))<1e-10_dp,'schur reconstruction')
  call check(maxval(abs(matmul(transpose(q),q)-ident))<1e-10_dp,'schur orthogonality')

  ar=reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp],[2,2])
  call schur_decomp(ar,tr,wrr,wir,qr,info)
  call check(info==0,'complex-pair schur')
  call check(maxval(abs(abs(wir)-1.0_dp))<1e-12_dp,'schur complex eigenvalues')

  x(:,1)=1.0_dp
  x(:,2)=[-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp,3.0_dp]
  y=2.0_dp+3.0_dp*x(:,2)
  call ols_fit_qr(x,y,rq,info)
  call check(info==0,'qr ols info')
  call check(maxval(abs(rq%coefficients-[2.0_dp,3.0_dp]))<1e-11_dp,'qr ols coefficients')
  call ols_fit_svd(x,y,rs,info)
  call check(info==0,'svd ols info')
  call check(maxval(abs(rs%coefficients-[2.0_dp,3.0_dp]))<1e-11_dp,'svd ols coefficients')

  asvd=x
  call svd_decomp(asvd,us,sv,vts,info)
  call check(info==0,'svd info')
  call check(maxval(abs(matmul(us,matmul(diag2(sv),vts))-asvd))<1e-10_dp,'svd reconstruction')

  utr=0.0_dp
  utr(1,1)=1.0_dp; utr(2,2)=2.0_dp; utr(3,3)=4.0_dp
  utr(1,2)=0.2_dp; utr(1,3)=-0.1_dp; utr(2,3)=0.3_dp
  call matrix_fun(utr,expfun,f,info)
  call check(info==0,'matrix fun info')
  call check(maxval(abs(matmul(f,utr)-matmul(utr,f)))<1e-10_dp,'matrix fun commutes')
  call check(maxval(abs([f(1,1)-exp(1.0_dp),f(2,2)-exp(2.0_dp),f(3,3)-exp(4.0_dp)]))<1e-12_dp,'matrix fun diagonal')

  xm(:,1)=[-1.0_dp,0.0_dp,1.0_dp,-1.0_dp,0.0_dp,1.0_dp,-0.5_dp,0.5_dp]
  xm(:,2)=[0.1_dp,1.2_dp,-0.7_dp,-1.1_dp,0.4_dp,0.8_dp,-0.3_dp,-0.4_dp]
  call mardia_test(xm,mr,info)
  call check(info==0,'mardia info')
  call check(mr%skewness>=0.0_dp .and. mr%kurtosis>0.0_dp,'mardia coefficients')
  call check(mr%skew_pvalue>=0.0_dp .and. mr%skew_pvalue<=1.0_dp,'mardia skew p')
  call check(mr%kurt_pvalue>=0.0_dp .and. mr%kurt_pvalue<=1.0_dp,'mardia kurt p')

  call harris_test(xm,'Wald',hr,info)
  call check(info==0,'harris wald info')
  call check(hr%statistic>=0.0_dp .and. hr%pvalue>=0.0_dp .and. hr%pvalue<=1.0_dp,'harris wald')
  call harris_test(xm,'log',hr,info)
  call check(info==0,'harris log info')
  call harris_test(xm,'robust',hr,info)
  call check(info==0,'harris robust info')
  call harris_test(xm,'log-robust',hr,info)
  call check(info==0,'harris log robust info')

  print '(a)', 'test_parity_v02: PASS'
contains
  function expfun(z) result(v)
    real(dp),intent(in)::z
    real(dp)::v
    v=exp(z)
  end function
  function diag2(s) result(d)
    real(dp),intent(in)::s(:)
    real(dp)::d(size(s),size(s))
    integer::j
    d=0.0_dp
    do j=1,size(s); d(j,j)=s(j); end do
  end function
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print '(a)', 'FAIL: '//trim(msg)
      error stop 1
    end if
  end subroutine
end program
