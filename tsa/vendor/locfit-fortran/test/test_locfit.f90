program test_locfit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use locfit
  implicit none
  integer :: failures
  failures=0
  call test_math(failures)
  call test_kernels(failures)
  call test_polynomial_fit(failures)
  call test_multivariate_fit(failures)
  call test_binomial_poisson(failures)
  call test_density(failures)
  call test_bandwidth(failures)
  call test_interpolation(failures)
  call test_diagnostics(failures)
  call test_auto_scale(failures)
  call test_robust(failures)
  if(failures/=0)then
    write(*,'(a,i0)') 'test_locfit: FAIL ',failures
    error stop 1
  end if
  write(*,'(a)') 'test_locfit: PASS'
contains
  subroutine check_close(name,got,want,tol,failures)
    character(*),intent(in)::name
    real(dp),intent(in)::got,want,tol
    integer,intent(inout)::failures
    if(.not.ieee_is_finite(got) .or. abs(got-want)>tol)then
      write(*,'(a,2(1x,es24.16),a,es12.4)') trim(name)//':',got,want,' tol=',tol
      failures=failures+1
    end if
  end subroutine check_close

  subroutine check_true(name,ok,failures)
    character(*),intent(in)::name
    logical,intent(in)::ok
    integer,intent(inout)::failures
    if(.not.ok)then
      write(*,'(a)') trim(name)//': false'
      failures=failures+1
    end if
  end subroutine check_true

  subroutine test_math(failures)
    integer,intent(inout)::failures
    call check_close('normal cdf zero',normal_cdf(0.0_dp),0.5_dp,2.0e-16_dp,failures)
    call check_close('gamma P(1,x)',gamma_p(0.7_dp,1.0_dp),1.0_dp-exp(-0.7_dp),2.0e-14_dp,failures)
    call check_close('beta I(.,1,1)',beta_i(0.37_dp,1.0_dp,1.0_dp),0.37_dp,2.0e-14_dp,failures)
    call check_close('expit large',expit(1000.0_dp),1.0_dp,0.0_dp,failures)
  end subroutine test_math

  subroutine test_kernels(failures)
    integer,intent(inout)::failures
    real(dp)::tf(19)
    integer::nc,st
    call check_close('tricube center',kernel_weight(0.0_dp,wtcub),1.0_dp,0.0_dp,failures)
    call check_close('tricube edge',kernel_weight(1.0_dp,wtcub),0.0_dp,0.0_dp,failures)
    call check_close('epan half',kernel_weight(0.5_dp,wepan),0.75_dp,1.0e-15_dp,failures)
    call check_close('gaussian integral',kernel_integral_moment(1,ker=wgaus), &
      sqrt(2.0_dp*pi)/gfact,2.0e-14_dp,failures)
    call check_close('rect convolution',kernel_convolution(0.5_dp,wrect),1.5_dp,1.0e-15_dp,failures)
    call kernel_taylor(0.0_dp,wtcub,tf,nc,st)
    call check_true('tricube taylor status',st==lf_ok .and. nc==10,failures)
    call check_close('tricube taylor z3',tf(4),-3.0_dp,0.0_dp,failures)
  end subroutine test_kernels

  subroutine test_polynomial_fit(failures)
    integer,intent(inout)::failures
    real(dp)::x(9,1),y(9),xe(3,1),der
    type(locfit_options)::op
    type(locfit_result)::fit
    integer::i,st
    do i=1,9
      x(i,1)=-2.0_dp+0.5_dp*real(i-1,dp)
      y(i)=1.0_dp+2.0_dp*x(i,1)+3.0_dp*x(i,1)**2
    end do
    xe(:,1)=[-1.3_dp,0.25_dp,1.1_dp]
    op%degree=2;op%kernel=wparm;op%nn=0.0_dp;op%h=1.0_dp
    op%family=tgaus;op%link=lident
    call locfit_fit(x,y,xe,fit,op)
    do i=1,3
      call check_true('quadratic status',fit%status(i)==lf_ok,failures)
      call check_close('quadratic exact fit',fit%fit(i), &
        1.0_dp+2.0_dp*xe(i,1)+3.0_dp*xe(i,1)**2,2.0e-11_dp,failures)
    end do
    call locfit_derivative_at(x,y,[0.25_dp],[1],der,st,op)
    call check_true('derivative status',st==lf_ok,failures)
    call check_close('quadratic derivative',der,3.5_dp,2.0e-11_dp,failures)
  end subroutine test_polynomial_fit

  subroutine test_multivariate_fit(failures)
    integer,intent(inout)::failures
    real(dp)::x(12,2),y(12),xe(2,2)
    type(locfit_options)::op
    type(locfit_result)::fit
    integer::i,j,k
    k=0
    do i=0,3
      do j=0,2
        k=k+1
        x(k,:)=[real(i,dp)-1.0_dp,real(j,dp)-1.0_dp]
        y(k)=2.0_dp+1.5_dp*x(k,1)-0.75_dp*x(k,2)+0.4_dp*x(k,1)**2+ &
          0.3_dp*x(k,1)*x(k,2)-0.2_dp*x(k,2)**2
      end do
    end do
    xe(1,:)=[0.2_dp,-0.4_dp];xe(2,:)=[1.1_dp,0.3_dp]
    op%degree=2;op%kernel=wparm;op%nn=0.0_dp;op%h=1.0_dp
    op%family=tgaus;op%link=lident;op%kernel_type=ksph
    call locfit_fit(x,y,xe,fit,op)
    do i=1,2
      call check_true('2d status',fit%status(i)==lf_ok,failures)
      call check_close('2d exact fit',fit%fit(i),2.0_dp+1.5_dp*xe(i,1)-0.75_dp*xe(i,2)+ &
        0.4_dp*xe(i,1)**2+0.3_dp*xe(i,1)*xe(i,2)-0.2_dp*xe(i,2)**2,5.0e-11_dp,failures)
    end do
  end subroutine test_multivariate_fit

  subroutine test_binomial_poisson(failures)
    integer,intent(inout)::failures
    real(dp)::x(5,1),yb(5),yp(5),xe(1,1)
    type(locfit_options)::op
    type(locfit_result)::fit
    integer::i
    do i=1,5;x(i,1)=real(i,dp);end do
    xe(1,1)=0.0_dp
    yb=[0.0_dp,1.0_dp,0.0_dp,1.0_dp,1.0_dp]
    op%degree=0;op%kernel=wparm;op%nn=0.0_dp;op%h=1.0_dp
    op%family=tlogt;op%link=llogit;op%compute_vcov=.false.
    call locfit_fit(x,yb,xe,fit,op)
    call check_true('binomial status',fit%status(1)==lf_ok,failures)
    call check_close('binomial intercept mean',fit%fit(1),0.6_dp,2.0e-11_dp,failures)
    yp=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    op%family=tpois;op%link=llog
    call locfit_fit(x,yp,xe,fit,op)
    call check_true('poisson status',fit%status(1)==lf_ok,failures)
    call check_close('poisson intercept mean',fit%fit(1),2.0_dp,2.0e-10_dp,failures)
  end subroutine test_binomial_poisson

  subroutine test_density(failures)
    integer,intent(inout)::failures
    real(dp)::x(7),pts(4),d1(4),d2(4),se(4)
    integer::st(4),i
    x=[-1.7_dp,-0.8_dp,-0.2_dp,0.0_dp,0.35_dp,0.9_dp,1.6_dp]
    pts=[-0.7_dp,0.0_dp,0.4_dp,1.2_dp]
    do i=1,4
      d1(i)=kernel_density_1d(x,pts(i),0.75_dp,wgaus)
    end do
    call local_density_1d(x,pts,0.75_dp,d2,degree=0,ker=wgaus,link=llog,status=st,se=se,nint=800)
    do i=1,4
      call check_true('density status',st(i)==lf_ok,failures)
      call check_close('degree0 density equals kde',d2(i),d1(i),2.0e-9_dp,failures)
      call check_true('density se finite',ieee_is_finite(se(i)) .and. se(i)>=0.0_dp,failures)
    end do
  end subroutine test_density

  subroutine test_bandwidth(failures)
    integer,intent(inout)::failures
    real(dp)::x(12),b(3),r(3)
    integer::m(3)
    x=[-2.0_dp,-1.5_dp,-1.1_dp,-0.7_dp,-0.2_dp,0.0_dp,0.2_dp,0.6_dp,0.9_dp,1.3_dp,1.7_dp,2.1_dp]
    call kde_criterion(x,0.7_dp,kde_lscv,wgaus,res=r)
    call check_true('lscv criterion finite',all(ieee_is_finite(r)),failures)
    m=[kde_lscv,kde_bcv,kde_gkk]
    call kde_select(x,0.25_dp,1.4_dp,m,wgaus,b)
    call check_true('bandwidths positive',all(ieee_is_finite(b)) .and. all(b>0.0_dp),failures)
  end subroutine test_bandwidth

  subroutine test_interpolation(failures)
    integer,intent(inout)::failures
    real(dp)::val(2,2),z
    val=0.0_dp;val(1,:)=[0.0_dp,0.0_dp];val(2,:)=[1.0_dp,2.0_dp]
    z=rectcell_interp([0.25_dp],val,[0.0_dp],[1.0_dp],2)
    call check_close('hermite x^2',z,0.0625_dp,2.0e-15_dp,failures)
    call check_close('cubic direct x^2',cubic_interp(0.25_dp,0.0_dp,1.0_dp,0.0_dp,2.0_dp), &
      0.0625_dp,2.0e-15_dp,failures)
  end subroutine test_interpolation

  subroutine test_diagnostics(failures)
    integer,intent(inout)::failures
    real(dp)::v,mrl(4)
    integer::st
    v=locfit_residual(3.0_dp,1.0_dp,2.0_dp,tgaus,lident,rraw,st)
    call check_true('raw residual status',st==lf_ok,failures)
    call check_close('raw residual',v,1.0_dp,2.0e-15_dp,failures)
    call check_close('aic formula',aic_score(-10.0_dp,3.0_dp),26.0_dp,0.0_dp,failures)
    call check_close('gcv formula',gcv_score(20,-10.0_dp,3.0_dp),400.0_dp/289.0_dp,2.0e-15_dp,failures)
    call check_close('cp formula',cp_score(20,-10.0_dp,3.0_dp,2.0_dp),-4.0_dp,2.0e-15_dp,failures)
    call km_mean_residual_life([1.0_dp,2.0_dp,3.0_dp,4.0_dp], &
      [.false.,.true.,.false.,.true.],mrl)
    call check_true('km mrl finite',all(ieee_is_finite(mrl)) .and. all(mrl>=0.0_dp),failures)
    call check_close('event mrl zero 1',mrl(1),0.0_dp,0.0_dp,failures)
    call check_close('event mrl zero 3',mrl(3),0.0_dp,0.0_dp,failures)
  end subroutine test_diagnostics

  subroutine test_auto_scale(failures)
    integer,intent(inout)::failures
    real(dp)::x(4,2),s(2)
    x(:,1)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    x(:,2)=[5.0_dp,5.0_dp,5.0_dp,5.0_dp]
    call automatic_scale(x,s)
    call check_close('automatic sd',s(1),sqrt(5.0_dp/3.0_dp),2.0e-15_dp,failures)
    call check_close('constant scale fallback',s(2),1.0_dp,0.0_dp,failures)
  end subroutine test_auto_scale

  subroutine test_robust(failures)
    integer,intent(inout)::failures
    real(dp)::x(7,1),y(7),xe(1,1)
    type(locfit_options)::op
    type(locfit_result)::ordinary,rob
    integer::i
    do i=1,7;x(i,1)=real(i,dp);end do
    y=[0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,100.0_dp]
    xe(1,1)=4.0_dp
    op%degree=0;op%kernel=wparm;op%nn=0.0_dp;op%h=1.0_dp
    op%family=tgaus;op%link=lident;op%compute_vcov=.false.
    call locfit_fit(x,y,xe,ordinary,op)
    call locfit_robust_fit(x,y,xe,rob,op,iterations=3)
    call check_true('robust status',rob%status(1)==lf_ok,failures)
    call check_true('robust outlier resistance',abs(rob%fit(1))<abs(ordinary%fit(1)),failures)
  end subroutine test_robust
end program test_locfit
