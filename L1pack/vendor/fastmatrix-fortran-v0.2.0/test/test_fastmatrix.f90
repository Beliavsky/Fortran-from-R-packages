program test_fastmatrix
  use fastmatrix
  implicit none
  real(dp) :: a(3,3), l(3,3), u(3,3), rec(3,3), inva(3,3)
  real(dp) :: eval(3), evec(3,3), x(3), b(3), sol(3), s(3,3), w(3,3)
  real(dp) :: xx(6,2), y(6), coef(2), sphere(500,3), ball(500,3), mu(3), sig(3,3), samp(100,3)
  real(dp) :: km(3,3), ctrl(3,2), tt(3), curve(3,2), dmat(3,3), dout(3,3)
  integer :: piv(3), info, it, i
  type(ols_result) :: fit
  type(moments_result) :: mm

  a = reshape([4.0_dp,1.0_dp,1.0_dp, 1.0_dp,3.0_dp,0.5_dp, 1.0_dp,0.5_dp,2.0_dp],[3,3])
  call lu_decomp(a,l,u,piv,info)
  if(info/=0) error stop 'LU failed'
  rec=matmul(l,u)
  if(maxval(abs(rec-a(piv,:)))>1e-10_dp) error stop 'LU identity failed'

  call lu_inverse(a,inva,info)
  if(maxval(abs(matmul(a,inva)-eye(3)))>1e-10_dp) error stop 'inverse failed'

  call jacobi_eigen(a,eval,evec,info=info)
  if(info/=0) error stop 'jacobi failed'
  do i=1,3
    if(maxval(abs(matmul(a,evec(:,i))-eval(i)*evec(:,i)))>1e-8_dp) error stop 'eigen residual failed'
  end do

  x=[1.0_dp,2.0_dp,3.0_dp]
  b=matmul(a,x)
  call cg_solve(a,b,sol,1e-12_dp,100,it)
  if(maxval(abs(sol-x))>1e-8_dp) error stop 'CG failed'
  call seidel_solve(a,b,sol,1e-12_dp,1000,it)
  if(maxval(abs(sol-x))>1e-8_dp) error stop 'Seidel failed'

  call matrix_sqrt(a,s,info)
  if(maxval(abs(matmul(s,s)-a))>1e-7_dp) error stop 'matrix sqrt failed'
  call whitening(a,w,info)
  if(maxval(abs(matmul(matmul(w,a),transpose(w))-eye(3)))>1e-7_dp) error stop 'whitening failed'

  if(maxval(abs(matmul(commutation(2,3),vec(reshape([1._dp,2._dp,3._dp,4._dp,5._dp,6._dp],[2,3]))) &
      -vec(transpose(reshape([1._dp,2._dp,3._dp,4._dp,5._dp,6._dp],[2,3])))))>1e-12_dp) error stop 'commutation failed'

  xx(:,1)=1.0_dp
  xx(:,2)=[0._dp,1._dp,2._dp,3._dp,4._dp,5._dp]
  coef=[2.0_dp,-0.5_dp]
  y=matmul(xx,coef)
  call ols_fit(xx,y,fit,info)
  if(maxval(abs(fit%coefficients-coef))>1e-10_dp) error stop 'OLS failed'
  call ridge_fit(xx,y,0.0_dp,sol(1:2))
  if(maxval(abs(sol(1:2)-coef))>1e-10_dp) error stop 'ridge failed'

  call rsphere(500,3,sphere)
  do i=1,500
    if(abs(sum(sphere(i,:)**2)-1.0_dp)>1e-10_dp) error stop 'sphere RNG failed'
  end do
  call rball(500,3,ball)
  do i=1,500
    if(sum(ball(i,:)**2)>1.0_dp+1e-12_dp) error stop 'ball RNG failed'
  end do

  mu=0.0_dp; sig=eye(3)
  call rmnorm(100,mu,sig,samp,info)
  if(info/=0) error stop 'rmnorm failed'

  mm=moments([1._dp,2._dp,3._dp,4._dp])
  if(abs(mm%mean-2.5_dp)>1e-12_dp) error stop 'moments failed'
  if(abs(ccc([1._dp,2._dp,3._dp],[1._dp,2._dp,3._dp])-1.0_dp)>1e-12_dp) error stop 'ccc failed'
  if(abs(pchi(qchi(0.8_dp,4.0_dp),4.0_dp)-0.8_dp)>1e-9_dp) error stop 'chi inversion failed'

  call krylov(a,[1._dp,0._dp,0._dp],3,km)
  if(maxval(abs(km(:,2)-matmul(a,km(:,1))))>1e-12_dp) error stop 'krylov failed'
  ctrl=reshape([0._dp,0._dp, 1._dp,2._dp, 2._dp,0._dp],[3,2])
  tt=[0._dp,0.5_dp,1._dp]
  call bezier_curve(ctrl,tt,curve)
  if(maxval(abs(curve(1,:)-ctrl(1,:)))>1e-12_dp .or. maxval(abs(curve(3,:)-ctrl(3,:)))>1e-12_dp) error stop 'bezier failed'

  dmat=reshape([0._dp,1._dp,10._dp,1._dp,0._dp,2._dp,10._dp,2._dp,0._dp],[3,3])
  call floyd_warshall(dmat,dout)
  if(abs(dout(1,3)-3._dp)>1e-12_dp) error stop 'floyd failed'

  print '(a)', 'test_fastmatrix: PASS'
end program
