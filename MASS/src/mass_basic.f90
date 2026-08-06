! SPDX-License-Identifier: GPL-3.0-only
module mass_basic
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rrcov_kinds, only : dp
  use rrcov_types, only : rrcov_success
  use rrcov_linalg, only : general_inverse, symmetric_eigen, identity_matrix
  use rrcov_random, only : seed_random, random_normal
  use mass_types, only : kde2d_result, mass_success, mass_invalid_argument, mass_dimension_error, mass_singular
  use mass_math, only : pi_dp, normal_pdf, sample_variance, type7_quantile, matrix_pseudoinverse
  implicit none
  private
  public :: ginv, null_space, mvrnorm, kde2d, bandwidth_nrd, nclass_freq
  public :: successive_difference_contrasts, adaptive_area, beta_kernel
  public :: rational_approximation, rational_values, write_matrix
  public :: con2tr

  abstract interface
    function scalar_function(x) result(value)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: value
    end function scalar_function
  end interface

contains

  function ginv(x, tolerance, status) result(value)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in), optional :: tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: value(:, :)
    value = matrix_pseudoinverse(x, tolerance, status)
  end function ginv

  subroutine null_space(m, basis, status, tolerance)
    real(dp), intent(in) :: m(:, :)
    real(dp), allocatable, intent(out) :: basis(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: gram(:, :), values(:), vectors(:, :)
    real(dp) :: tol, vmax
    integer :: st, k
    if (size(m,2) == 0) then
      allocate(basis(0,0)); status=mass_invalid_argument; return
    end if
    gram = matmul(transpose(m),m)
    call symmetric_eigen(gram,values,vectors,st)
    vmax=max(1.0_dp,maxval(abs(values)))
    tol=sqrt(epsilon(1.0_dp)); if(present(tolerance)) tol=tolerance
    k=count(values <= tol*vmax)
    allocate(basis(size(m,2),k))
    if(k>0) basis=vectors(:,size(values)-k+1:size(values))
    status=merge(mass_success,mass_singular,st==rrcov_success)
  end subroutine null_space

  subroutine mvrnorm(n, mu, sigma, samples, status, tolerance, empirical, seed)
    integer, intent(in) :: n
    real(dp), intent(in) :: mu(:), sigma(:, :)
    real(dp), allocatable, intent(out) :: samples(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: tolerance
    logical, intent(in), optional :: empirical
    integer, intent(in), optional :: seed
    real(dp), allocatable :: values(:),vectors(:,:),z(:,:),zc(:,:),covz(:,:),rootinv(:,:)
    real(dp) :: tol
    logical :: emp
    integer :: p,i,j,st
    p=size(mu); tol=1.0e-6_dp; if(present(tolerance)) tol=tolerance
    emp=.false.; if(present(empirical)) emp=empirical
    if(n<1 .or. size(sigma,1)/=p .or. size(sigma,2)/=p) then
      allocate(samples(0,0)); status=mass_dimension_error; return
    end if
    call symmetric_eigen(sigma,values,vectors,st)
    if(minval(values) < -tol*max(1.0_dp,abs(values(1)))) then
      allocate(samples(0,0)); status=mass_invalid_argument; return
    end if
    call seed_random(seed)
    allocate(z(n,p))
    do i=1,n; do j=1,p; z(i,j)=random_normal(); end do; end do
    if(emp .and. n>p) then
      zc=z-spread(sum(z,dim=1)/real(n,dp),1,n)
      covz=matmul(transpose(zc),zc)/real(n-1,dp)
      call inverse_symmetric_sqrt(covz,rootinv,st)
      z=matmul(zc,rootinv)
    end if
    allocate(samples(n,p))
    samples=matmul(z,matmul(diag_sqrt(max(values,0.0_dp)),transpose(vectors)))
    samples=samples+spread(mu,1,n)
    status=mass_success
  contains
    function diag_sqrt(v) result(d)
      real(dp),intent(in)::v(:); real(dp)::d(size(v),size(v)); integer::q
      d=0.0_dp; do q=1,size(v); d(q,q)=sqrt(v(q)); end do
    end function diag_sqrt
    subroutine inverse_symmetric_sqrt(a,b,ist)
      real(dp),intent(in)::a(:,:); real(dp),allocatable,intent(out)::b(:,:); integer,intent(out)::ist
      real(dp),allocatable::ev(:),vec(:,:),d(:,:); integer::q
      call symmetric_eigen(a,ev,vec,ist); allocate(d(size(ev),size(ev))); d=0.0_dp
      do q=1,size(ev); d(q,q)=1.0_dp/sqrt(max(ev(q),1.0e-12_dp)); end do
      b=matmul(vec,matmul(d,transpose(vec)))
    end subroutine inverse_symmetric_sqrt
  end subroutine mvrnorm

  subroutine kde2d(x,y,result,h,n_grid,limits)
    real(dp),intent(in)::x(:),y(:)
    type(kde2d_result),intent(out)::result
    real(dp),intent(in),optional::h(:),limits(:)
    integer,intent(in),optional::n_grid(:)
    real(dp)::bw(2),lims(4),dx,dy
    integer::ng(2),i,j,k,n
    n=size(x)
    if(n==0 .or. size(y)/=n .or. any(.not.ieee_is_finite(x)) .or. any(.not.ieee_is_finite(y))) then
      result%status=mass_invalid_argument; return
    end if
    ng=[25,25]; if(present(n_grid)) then
      if(size(n_grid)==1) ng=n_grid(1); if(size(n_grid)>=2) ng=n_grid(1:2)
    end if
    lims=[minval(x),maxval(x),minval(y),maxval(y)]; if(present(limits)) then
      if(size(limits)>=4) lims=limits(1:4)
    end if
    bw=[bandwidth_nrd(x),bandwidth_nrd(y)]; if(present(h)) then
      if(size(h)==1) bw=h(1); if(size(h)>=2) bw=h(1:2)
    end if
    if(any(bw<=0.0_dp) .or. any(ng<2)) then; result%status=mass_invalid_argument; return; end if
    bw=bw/4.0_dp
    allocate(result%x_grid(ng(1)),result%y_grid(ng(2)),result%density(ng(1),ng(2)))
    dx=(lims(2)-lims(1))/real(ng(1)-1,dp); dy=(lims(4)-lims(3))/real(ng(2)-1,dp)
    result%x_grid=[(lims(1)+real(i-1,dp)*dx,i=1,ng(1))]
    result%y_grid=[(lims(3)+real(j-1,dp)*dy,j=1,ng(2))]
    result%density=0.0_dp
    do j=1,ng(2); do i=1,ng(1)
      do k=1,n
        result%density(i,j)=result%density(i,j)+normal_pdf((result%x_grid(i)-x(k))/bw(1))* &
          normal_pdf((result%y_grid(j)-y(k))/bw(2))
      end do
      result%density(i,j)=result%density(i,j)/(real(n,dp)*bw(1)*bw(2))
    end do; end do
    result%status=mass_success
  end subroutine kde2d

  function bandwidth_nrd(x) result(h)
    real(dp),intent(in)::x(:)
    real(dp)::h,iqr,s
    if(size(x)<2) then; h=0.0_dp; return; end if
    iqr=(type7_quantile(x,0.75_dp)-type7_quantile(x,0.25_dp))/1.34_dp
    s=sqrt(max(0.0_dp,sample_variance(x)))
    h=4.0_dp*1.06_dp*min(s,iqr)*real(size(x),dp)**(-0.2_dp)
    if(h<=0.0_dp) h=4.0_dp*1.06_dp*max(s,iqr)*real(size(x),dp)**(-0.2_dp)
  end function bandwidth_nrd

  function nclass_freq(x) result(number)
    real(dp),intent(in)::x(:)
    integer::number
    real(dp)::h
    if(size(x)<2) then; number=1; return; end if
    h=2.15_dp*sqrt(max(0.0_dp,sample_variance(x)))*real(size(x),dp)**(-0.2_dp)
    if(h<=0.0_dp) then; number=1; else; number=ceiling((maxval(x)-minval(x))/h); end if
    number=max(1,number)
  end function nclass_freq

  function successive_difference_contrasts(n, contrasts) result(value)
    integer,intent(in)::n
    logical,intent(in),optional::contrasts
    real(dp),allocatable::value(:,:)
    logical::usec
    integer::i,j
    usec=.true.; if(present(contrasts)) usec=contrasts
    if(n<2) then; allocate(value(0,0)); return; end if
    if(.not.usec) then; value=identity_matrix(n); return; end if
    allocate(value(n,n-1))
    do j=1,n-1
      do i=1,n
        value(i,j)=real(j,dp)/real(n,dp)
        if(i<=j) value(i,j)=value(i,j)-1.0_dp
      end do
    end do
  end function successive_difference_contrasts

  recursive function adaptive_area(fun,a,b,limit,eps,status) result(value)
    procedure(scalar_function)::fun
    real(dp),intent(in)::a,b
    integer,intent(in),optional::limit
    real(dp),intent(in),optional::eps
    integer,intent(out),optional::status
    real(dp)::value,d,fa,fb,fd,a1,a2,tol
    integer::lim,st1,st2
    lim=10; if(present(limit)) lim=limit
    tol=1.0e-5_dp; if(present(eps)) tol=eps
    fa=fun(a); fb=fun(b); d=0.5_dp*(a+b); fd=fun(d)
    a1=0.5_dp*(fa+fb)*(b-a); a2=(fa+4.0_dp*fd+fb)*(b-a)/6.0_dp
    if(abs(a1-a2)<tol .or. lim==0) then
      value=a2; if(present(status)) status=merge(mass_success,mass_no_convergence_local(),lim>0); return
    end if
    value=adaptive_area(fun,a,d,lim-1,tol,st1)+adaptive_area(fun,d,b,lim-1,tol,st2)
    if(present(status)) status=max(st1,st2)
  contains
    pure integer function mass_no_convergence_local()
      mass_no_convergence_local=4
    end function mass_no_convergence_local
  end function adaptive_area

  pure elemental function beta_kernel(x,alpha,beta) result(value)
    real(dp),intent(in)::x,alpha,beta
    real(dp)::value
    if(x<0.0_dp .or. x>1.0_dp) then; value=0.0_dp
    else; value=x**(alpha-1.0_dp)*(1.0_dp-x)**(beta-1.0_dp); end if
  end function beta_kernel

  subroutine rational_approximation(x,numerator,denominator,cycles,max_denominator)
    real(dp),intent(in)::x(:)
    integer(kind=selected_int_kind(18)),allocatable,intent(out)::numerator(:),denominator(:)
    integer,intent(in),optional::cycles,max_denominator
    integer(kind=selected_int_kind(18))::p0,p1,p2,q0,q1,q2,a
    real(dp)::z,frac
    integer::i,k,ncy,mden
    ncy=10; if(present(cycles)) ncy=cycles
    mden=2000; if(present(max_denominator)) mden=max_denominator
    allocate(numerator(size(x)),denominator(size(x)))
    do i=1,size(x)
      if(.not.ieee_is_finite(x(i))) then; numerator(i)=0; denominator(i)=0; cycle; end if
      z=abs(x(i)); p0=0; q0=1; p1=1; q1=0
      do k=1,ncy+1
        a=int(floor(z),kind(a)); p2=a*p1+p0; q2=a*q1+q0
        if(q2>mden .or. q2<=0) exit
        p0=p1; q0=q1; p1=p2; q1=q2
        frac=z-real(a,dp); if(frac<=1.0_dp/real(mden,dp)) exit
        z=1.0_dp/frac
      end do
      numerator(i)=merge(-p1,p1,x(i)<0.0_dp); denominator(i)=q1
    end do
  end subroutine rational_approximation

  function rational_values(x,cycles,max_denominator) result(value)
    real(dp),intent(in)::x(:)
    integer,intent(in),optional::cycles,max_denominator
    real(dp)::value(size(x))
    integer(kind=selected_int_kind(18)),allocatable::num(:),den(:)
    integer::i
    call rational_approximation(x,num,den,cycles,max_denominator)
    do i=1,size(x)
      if(den(i)==0) then; value(i)=x(i); else; value(i)=real(num(i),dp)/real(den(i),dp); end if
    end do
  end function rational_values

  subroutine write_matrix(x,filename,separator,status)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in)::filename
    character(len=*),intent(in),optional::separator
    integer,intent(out),optional::status
    character(len=16)::sep
    integer::u,ios,i,j
    sep=" "; if(present(separator)) sep=separator
    open(newunit=u,file=filename,status='replace',action='write',iostat=ios)
    if(ios/=0) then; if(present(status)) status=mass_invalid_argument; return; end if
    do i=1,size(x,1)
      do j=1,size(x,2)
        if(j>1) write(u,'(a)',advance='no') trim(sep)
        write(u,'(es24.16)',advance='no') x(i,j)
      end do
      write(u,*)
    end do
    close(u); if(present(status)) status=mass_success
  end subroutine write_matrix

  subroutine con2tr(x,y,z,table,status)
    real(dp),intent(in)::x(:),y(:),z(:,:)
    real(dp),allocatable,intent(out)::table(:,:)
    integer,intent(out)::status
    integer::i,j,k
    if(size(z,1)/=size(x) .or. size(z,2)/=size(y)) then
      allocate(table(0,0)); status=mass_dimension_error; return
    end if
    allocate(table(size(x)*size(y),3)); k=0
    do j=1,size(y); do i=1,size(x); k=k+1; table(k,:)=[x(i),y(j),z(i,j)]; end do; end do
    status=mass_success
  end subroutine con2tr

end module mass_basic
