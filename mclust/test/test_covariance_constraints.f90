program test_covariance_constraints
  use mclust
  use mclust_linalg, only : symmetric_eigen
  implicit none
  integer,parameter::n=150,d=3,g=3
  real(dp)::x(n,d),tol
  real(dp),allocatable::z0(:,:),ev(:,:),v(:,:),evals(:)
  type(mclust_fit)::fit
  character(len=3),parameter::mods(14)=[character(len=3) :: &
    'EII','VII','EEI','VEI','EVI','VVI','EEE', &
    'EVE','VEE','VVE','EEV','VEV','EVV','VVV']
  integer::i,k,j,st
  real(dp)::det(g),scale(g),base(d),err

  do i=1,n/3
    x(i,:)=[-3.0_dp+0.35_dp*sin(0.31_dp*i), -1.5_dp+0.18_dp*cos(0.27_dp*i), 0.5_dp+0.12_dp*sin(0.51_dp*i)]
  end do
  do i=n/3+1,2*n/3
    x(i,:)=[0.5_dp+0.15_dp*sin(0.37_dp*i), 3.0_dp+0.45_dp*cos(0.19_dp*i), -2.0_dp+0.22_dp*sin(0.47_dp*i)]
  end do
  do i=2*n/3+1,n
    x(i,:)=[4.0_dp+0.28_dp*sin(0.29_dp*i), -2.5_dp+0.30_dp*cos(0.33_dp*i), 3.0_dp+0.40_dp*sin(0.23_dp*i)]
  end do
  call hc_responsibilities(x,g,z0,'VVV',st); if(st/=0)error stop 'hc'
  tol=2e-5_dp
  do j=1,size(mods)
    call fit_model(x,g,mods(j),fit,z_init=z0)
    if(fit%status<0) then; print *,'failed ',mods(j),fit%status; error stop 'fit'; end if
    do k=1,g; det(k)=det3(fit%sigma(:,:,k)); end do
    select case(mods(j))
    case('EII')
      if(maxval(abs(fit%sigma(:,:,1)-diag3(fit%sigma(1,1,1))))>tol) error stop 'EII spherical'
      if(maxval(abs(fit%sigma-spread(fit%sigma(:,:,1),3,g)))>tol) error stop 'EII equal'
    case('VII')
      do k=1,g
        if(maxval(abs(fit%sigma(:,:,k)-diag3(fit%sigma(1,1,k))))>tol)error stop 'VII spherical'
      end do
    case('EEI')
      do k=1,g
        if(maxval(abs(fit%sigma(:,:,k)-diag3v([fit%sigma(1,1,k),fit%sigma(2,2,k),fit%sigma(3,3,k)])))>tol)error stop 'EEI diagonal'
      end do
      if(maxval(abs(fit%sigma-spread(fit%sigma(:,:,1),3,g)))>tol)error stop 'EEI equal'
    case('VEI')
      do k=1,g
        if(maxval(abs(fit%sigma(:,:,k)-diag3v([fit%sigma(1,1,k),fit%sigma(2,2,k),fit%sigma(3,3,k)])))>tol)error stop 'VEI diagonal'
        scale(k)=det(k)**(1.0_dp/3.0_dp)
      end do
      base=[fit%sigma(1,1,1),fit%sigma(2,2,1),fit%sigma(3,3,1)]/scale(1)
      do k=2,g
        err=maxval(abs([fit%sigma(1,1,k),fit%sigma(2,2,k),fit%sigma(3,3,k)]/scale(k)-base))
        if(err>tol)error stop 'VEI common shape'
      end do
    case('EVI')
      if(maxval(abs(det-det(1)))>tol*max(1.0_dp,abs(det(1))))error stop 'EVI equal volume'
      do k=1,g
        if(maxval(abs(fit%sigma(:,:,k)-diag3v([fit%sigma(1,1,k),fit%sigma(2,2,k),fit%sigma(3,3,k)])))>tol)error stop 'EVI diagonal'
      end do
    case('VVI')
      do k=1,g
        if(maxval(abs(fit%sigma(:,:,k)-diag3v([fit%sigma(1,1,k),fit%sigma(2,2,k),fit%sigma(3,3,k)])))>tol)error stop 'VVI diagonal'
      end do
    case('EEE')
      if(maxval(abs(fit%sigma-spread(fit%sigma(:,:,1),3,g)))>tol)error stop 'EEE equal'
    case('EVE')
      if(maxval(abs(det-det(1)))>tol*max(1.0_dp,abs(det(1))))error stop 'EVE equal volume'
      do k=2,g
        if(norm2(matmul(fit%sigma(:,:,1),fit%sigma(:,:,k)) - &
          matmul(fit%sigma(:,:,k),fit%sigma(:,:,1)))>1e-4_dp) &
          error stop 'EVE common orientation'
      end do
    case('VEE')
      do k=1,g; scale(k)=det(k)**(1.0_dp/3.0_dp); end do
      do k=2,g
        if(maxval(abs(fit%sigma(:,:,k)/scale(k)-fit%sigma(:,:,1)/scale(1)))>tol)error stop 'VEE common shape/orientation'
      end do
    case('VVE')
      do k=2,g
        if(norm2(matmul(fit%sigma(:,:,1),fit%sigma(:,:,k)) - &
          matmul(fit%sigma(:,:,k),fit%sigma(:,:,1)))>1e-4_dp) &
          error stop 'VVE common orientation'
      end do
    case('EEV')
      allocate(ev(d,g))
      do k=1,g
        call symmetric_eigen(fit%sigma(:,:,k),evals,v,st); if(st/=0)error stop 'eig'; ev(:,k)=evals
      end do
      do k=2,g; if(maxval(abs(ev(:,k)-ev(:,1)))>tol)error stop 'EEV common eigenvalues'; end do
      deallocate(ev)
    case('VEV')
      allocate(ev(d,g))
      do k=1,g
        call symmetric_eigen(fit%sigma(:,:,k),evals,v,st)
        if(st/=0)error stop 'eig'
        ev(:,k)=evals/(product(evals)**(1.0_dp/3.0_dp))
      end do
      do k=2,g; if(maxval(abs(ev(:,k)-ev(:,1)))>tol)error stop 'VEV common shape'; end do
      deallocate(ev)
    case('EVV')
      if(maxval(abs(det-det(1)))>tol*max(1.0_dp,abs(det(1))))error stop 'EVV equal volume'
    case('VVV')
      continue
    end select
  end do
  print *, 'test_covariance_constraints PASS'
contains
  pure real(dp) function det3(a)
    real(dp),intent(in)::a(3,3)
    det3=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  end function
  pure function diag3(a) result(m)
    real(dp),intent(in)::a; real(dp)::m(3,3); integer::ii
    m=0.0_dp; do ii=1,3;m(ii,ii)=a;end do
  end function
  pure function diag3v(a) result(m)
    real(dp),intent(in)::a(3); real(dp)::m(3,3); integer::ii
    m=0.0_dp; do ii=1,3;m(ii,ii)=a(ii);end do
  end function
end program test_covariance_constraints
