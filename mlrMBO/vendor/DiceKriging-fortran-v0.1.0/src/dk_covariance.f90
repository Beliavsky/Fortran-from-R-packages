! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! This translation is distributed under the same license choice; see
! LICENSE-GPL-2 and LICENSE-GPL-3 in the project root.
module dk_covariance
  use dk_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: cov_gauss = 1
  integer, parameter, public :: cov_exp = 2
  integer, parameter, public :: cov_matern32 = 3
  integer, parameter, public :: cov_matern52 = 4
  integer, parameter, public :: cov_powexp = 5

  type, public :: scaling_axis
    real(dp), allocatable :: knots(:)
    real(dp), allocatable :: eta(:)
  end type scaling_axis

  type, public :: covariance_model
    integer :: kind = cov_matern52
    logical :: iso = .false.
    logical :: scaling = .false.
    real(dp), allocatable :: range(:)
    real(dp), allocatable :: shape(:)
    real(dp) :: sd2 = 1.0_dp
    logical :: nugget_flag = .false.
    logical :: nugget_estim = .false.
    real(dp) :: nugget = 0.0_dp
    type(scaling_axis), allocatable :: axis(:)
  end type covariance_model

  public :: covariance_kind, covariance_name, covariance_param_count
  public :: covariance_matrix, covariance_cross, covariance_derivative
  public :: covariance_vector_dx, covariance_bounds
  public :: scaling_fun1d, scaling_fun1d_dx, scaling_grad1d, scaling_apply
  public :: get_cov_params, set_cov_params

contains

  integer function covariance_kind(name) result(k)
    character(len=*), intent(in) :: name
    select case (trim(adjustl(name)))
    case ('gauss'); k = cov_gauss
    case ('exp'); k = cov_exp
    case ('matern3_2'); k = cov_matern32
    case ('matern5_2'); k = cov_matern52
    case ('powexp'); k = cov_powexp
    case default; k = 0
    end select
  end function covariance_kind

  function covariance_name(k) result(name)
    integer, intent(in) :: k
    character(len=16) :: name
    select case (k)
    case (cov_gauss); name = 'gauss'
    case (cov_exp); name = 'exp'
    case (cov_matern32); name = 'matern3_2'
    case (cov_matern52); name = 'matern5_2'
    case (cov_powexp); name = 'powexp'
    case default; name = 'unknown'
    end select
  end function covariance_name

  integer function covariance_param_count(cov, d) result(n)
    type(covariance_model), intent(in) :: cov
    integer, intent(in) :: d
    integer :: j
    if (cov%scaling) then
      n = 0
      if (allocated(cov%axis)) then
        do j = 1, size(cov%axis)
          if (allocated(cov%axis(j)%eta)) n = n + size(cov%axis(j)%eta)
        end do
      end if
    else if (cov%iso) then
      n = 1
    else
      n = d
      if (cov%kind == cov_powexp) n = 2*d
    end if
  end function covariance_param_count

  subroutine get_cov_params(cov, d, p)
    type(covariance_model), intent(in) :: cov
    integer, intent(in) :: d
    real(dp), allocatable, intent(out) :: p(:)
    integer :: n, j, pos
    n = covariance_param_count(cov,d)
    allocate(p(n)); p = 0.0_dp
    if (cov%scaling) then
      pos = 1
      do j = 1, size(cov%axis)
        p(pos:pos+size(cov%axis(j)%eta)-1) = cov%axis(j)%eta
        pos = pos + size(cov%axis(j)%eta)
      end do
    else if (cov%iso) then
      p(1) = cov%range(1)
    else
      p(1:d) = cov%range(1:d)
      if (cov%kind == cov_powexp) p(d+1:2*d) = cov%shape(1:d)
    end if
  end subroutine get_cov_params

  subroutine set_cov_params(cov, d, p)
    type(covariance_model), intent(inout) :: cov
    integer, intent(in) :: d
    real(dp), intent(in) :: p(:)
    integer :: j, pos
    if (cov%scaling) then
      pos = 1
      do j = 1, size(cov%axis)
        cov%axis(j)%eta = p(pos:pos+size(cov%axis(j)%eta)-1)
        pos = pos + size(cov%axis(j)%eta)
      end do
    else if (cov%iso) then
      if (.not. allocated(cov%range)) allocate(cov%range(d))
      cov%range = p(1)
    else
      if (.not. allocated(cov%range)) allocate(cov%range(d))
      cov%range(1:d) = p(1:d)
      if (cov%kind == cov_powexp) then
        if (.not. allocated(cov%shape)) allocate(cov%shape(d))
        cov%shape(1:d) = p(d+1:2*d)
      end if
    end if
  end subroutine set_cov_params

  pure real(dp) function scaling_factor(kind) result(s)
    integer, intent(in) :: kind
    select case (kind)
    case (cov_gauss); s = sqrt(2.0_dp)/2.0_dp
    case (cov_matern32); s = sqrt(3.0_dp)
    case (cov_matern52); s = sqrt(5.0_dp)
    case default; s = 1.0_dp
    end select
  end function scaling_factor

  pure real(dp) function cov_pair_raw(kind, x1, x2, range, shape, var) result(v)
    integer, intent(in) :: kind
    real(dp), intent(in) :: x1(:), x2(:), range(:), var
    real(dp), intent(in), optional :: shape(:)
    integer :: j
    real(dp) :: s, e, scf
    scf = scaling_factor(kind)
    s = 0.0_dp
    select case (kind)
    case (cov_gauss)
      do j = 1, size(x1)
        e = (x1(j)-x2(j))/(range(j)/scf)
        s = s + e*e
      end do
      v = var*exp(-s)
    case (cov_exp)
      do j = 1, size(x1)
        s = s + abs(x1(j)-x2(j))/range(j)
      end do
      v = var*exp(-s)
    case (cov_matern32)
      do j = 1, size(x1)
        e = abs(x1(j)-x2(j))/(range(j)/scf)
        s = s + e - log(1.0_dp+e)
      end do
      v = var*exp(-s)
    case (cov_matern52)
      do j = 1, size(x1)
        e = abs(x1(j)-x2(j))/(range(j)/scf)
        s = s + e - log(1.0_dp+e+e*e/3.0_dp)
      end do
      v = var*exp(-s)
    case (cov_powexp)
      do j = 1, size(x1)
        s = s + (abs(x1(j)-x2(j))/range(j))**shape(j)
      end do
      v = var*exp(-s)
    case default
      v = 0.0_dp
    end select
  end function cov_pair_raw

  subroutine effective_ranges(cov, d, range, shape)
    type(covariance_model), intent(in) :: cov
    integer, intent(in) :: d
    real(dp), allocatable, intent(out) :: range(:), shape(:)
    allocate(range(d)); allocate(shape(d))
    if (cov%iso) then
      range = cov%range(1)
    else
      range = cov%range(1:d)
    end if
    shape = 1.0_dp
    if (cov%kind == cov_powexp .and. allocated(cov%shape)) shape = cov%shape(1:d)
  end subroutine effective_ranges

  subroutine covariance_matrix(cov, x, c, noise_var, include_nugget)
    type(covariance_model), intent(in) :: cov
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: c(:,:)
    real(dp), intent(in), optional :: noise_var(:)
    logical, intent(in), optional :: include_nugget
    real(dp), allocatable :: xx(:,:), range(:), shape(:)
    integer :: n, d, i, j
    logical :: add_nug
    n = size(x,1); d = size(x,2)
    if (cov%scaling) then
      call scaling_apply(x, cov%axis, xx)
    else
      xx = x
    end if
    call effective_ranges(cov,d,range,shape)
    allocate(c(n,n)); c = 0.0_dp
    do i = 1, n
      c(i,i) = cov%sd2
      do j = 1, i-1
        c(i,j) = cov_pair_raw(cov%kind,xx(i,:),xx(j,:),range,shape,cov%sd2)
        c(j,i) = c(i,j)
      end do
    end do
    add_nug = cov%nugget_flag
    if (present(include_nugget)) add_nug = include_nugget .and. cov%nugget_flag
    if (add_nug) then
      do i = 1, n
        c(i,i) = c(i,i) + cov%nugget
      end do
    else if (present(noise_var)) then
      do i = 1, n
        c(i,i) = c(i,i) + noise_var(i)
      end do
    end if
  end subroutine covariance_matrix

  subroutine covariance_cross(cov, x1, x2, c, include_nugget)
    type(covariance_model), intent(in) :: cov
    real(dp), intent(in) :: x1(:,:), x2(:,:)
    real(dp), allocatable, intent(out) :: c(:,:)
    logical, intent(in), optional :: include_nugget
    real(dp), allocatable :: a(:,:), b(:,:), range(:), shape(:)
    integer :: n1, n2, d, i, j
    logical :: add_nug
    n1=size(x1,1); n2=size(x2,1); d=size(x1,2)
    if (cov%scaling) then
      call scaling_apply(x1,cov%axis,a); call scaling_apply(x2,cov%axis,b)
    else
      a=x1; b=x2
    end if
    call effective_ranges(cov,d,range,shape)
    allocate(c(n1,n2))
    do j=1,n2
      do i=1,n1
        c(i,j)=cov_pair_raw(cov%kind,a(i,:),b(j,:),range,shape,cov%sd2)
      end do
    end do
    add_nug=.false.; if (present(include_nugget)) add_nug=include_nugget
    if (add_nug .and. cov%nugget_flag) then
      do j=1,n2
        do i=1,n1
          if (sum(abs(x1(i,:)-x2(j,:))) < 1.0e-15_dp) c(i,j)=c(i,j)+cov%nugget
        end do
      end do
    end if
  end subroutine covariance_cross

  subroutine covariance_derivative(cov, x, c0, k, dc)
    type(covariance_model), intent(in) :: cov
    real(dp), intent(in) :: x(:,:), c0(:,:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: dc(:,:)
    type(covariance_model) :: cp, cm
    real(dp), allocatable :: p(:), pp(:), pm(:), cplus(:,:), cminus(:,:)
    real(dp) :: h
    integer :: n, d, i, j, kk
    real(dp), allocatable :: xx(:,:), range(:), shape(:)
    real(dp) :: e, scf, dln, r
    n=size(x,1); d=size(x,2)
    allocate(dc(n,n)); dc=0.0_dp
    if (cov%scaling) then
      call get_cov_params(cov,d,p)
      if (k <= size(p)) then
        h=max(1.0e-7_dp,1.0e-5_dp*abs(p(k)))
        pp=p; pm=p; pp(k)=p(k)+h; pm(k)=max(1.0e-12_dp,p(k)-h)
        cp=cov; cm=cov
        call set_cov_params(cp,d,pp); call set_cov_params(cm,d,pm)
        call covariance_matrix(cp,x,cplus,include_nugget=.false.)
        call covariance_matrix(cm,x,cminus,include_nugget=.false.)
        dc=(cplus-cminus)/(pp(k)-pm(k))
      else if (k == size(p)+1) then
        dc=c0/cov%sd2
      end if
      return
    end if
    if (k == covariance_param_count(cov,d)+1) then
      dc=c0/cov%sd2
      return
    end if
    xx=x; call effective_ranges(cov,d,range,shape); scf=scaling_factor(cov%kind)
    do i=1,n
      do j=1,i-1
        if (cov%iso) then
          dln=0.0_dp
          do kk=1,d
            r=range(kk)
            select case(cov%kind)
            case(cov_gauss)
              e=(xx(i,kk)-xx(j,kk))/(r/scf); dln=dln+2.0_dp*e*e/r
            case(cov_exp)
              e=abs(xx(i,kk)-xx(j,kk))/r; dln=dln+e/r
            case(cov_matern32)
              e=abs(xx(i,kk)-xx(j,kk))/(r/scf); dln=dln+(e*e/(1.0_dp+e))/r
            case(cov_matern52)
              e=abs(xx(i,kk)-xx(j,kk))/(r/scf)
              dln=dln+((1.0_dp+e)*(e*e/3.0_dp)/(1.0_dp+e+e*e/3.0_dp))/r
            end select
          end do
        else
          kk=mod(k-1,d)+1; r=range(kk)
          e=abs(xx(i,kk)-xx(j,kk))/r
          select case(cov%kind)
          case(cov_gauss)
            e=(xx(i,kk)-xx(j,kk))/(r/scf); dln=2.0_dp*e*e/r
          case(cov_exp)
            dln=e/r
          case(cov_matern32)
            e=abs(xx(i,kk)-xx(j,kk))/(r/scf); dln=(e*e/(1.0_dp+e))/r
          case(cov_matern52)
            e=abs(xx(i,kk)-xx(j,kk))/(r/scf)
            dln=((1.0_dp+e)*(e*e/3.0_dp)/(1.0_dp+e+e*e/3.0_dp))/r
          case(cov_powexp)
            if (k <= d) then
              dln=e**shape(kk)*shape(kk)/r
            else
              if (abs(e) <= tiny(1.0_dp)) then
                dln=0.0_dp
              else
                dln=-(e**shape(kk))*log(e)
              end if
            end if
          end select
        end if
        dc(i,j)=dln*c0(i,j); dc(j,i)=dc(i,j)
      end do
    end do
  end subroutine covariance_derivative

  subroutine covariance_vector_dx(cov, x, design, c, grad)
    type(covariance_model), intent(in) :: cov
    real(dp), intent(in) :: x(:), design(:,:), c(:)
    real(dp), allocatable, intent(out) :: grad(:,:)
    real(dp), allocatable :: xx(:), xd(:,:), range(:), shape(:), gx(:)
    integer :: n,d,i,k
    real(dp) :: e,sgn,scf,u,v
    n=size(design,1); d=size(design,2)
    if (cov%scaling) then
      allocate(xx(d),gx(d))
      do k=1,d
        xx(k)=scaling_fun1d(x(k),cov%axis(k)%knots,cov%axis(k)%eta)
        gx(k)=scaling_fun1d_dx(x(k),cov%axis(k)%knots,cov%axis(k)%eta)
      end do
      call scaling_apply(design,cov%axis,xd)
    else
      allocate(xx(d), xd(n,d), gx(d))
      xx = x
      xd = design
      gx = 1.0_dp
    end if
    call effective_ranges(cov,d,range,shape); scf=scaling_factor(cov%kind)
    allocate(grad(n,d)); grad=0.0_dp
    do i=1,n
      do k=1,d
        e=xx(k)-xd(i,k)
        if (abs(e) <= tiny(1.0_dp)) cycle
        sgn=merge(1.0_dp,-1.0_dp,e>0.0_dp)
        select case(cov%kind)
        case(cov_gauss)
          grad(i,k)=c(i)*(-2.0_dp*e/(range(k)/scf)**2)
        case(cov_exp)
          grad(i,k)=c(i)*(-sgn/range(k))
        case(cov_matern32)
          u=abs(e)/(range(k)/scf)
          grad(i,k)=c(i)*(-sgn*u/(1.0_dp+u)/(range(k)/scf))
        case(cov_matern52)
          u=abs(e)/(range(k)/scf); v=u/3.0_dp
          grad(i,k)=c(i)*(-sgn*(1.0_dp+u)*v/(1.0_dp+u+u*v)/(range(k)/scf))
        case(cov_powexp)
          u=abs(e)/range(k)
          grad(i,k)=c(i)*(-sgn*u**(shape(k)-1.0_dp)*shape(k)/range(k))
        end select
        grad(i,k)=grad(i,k)*gx(k)
      end do
    end do
  end subroutine covariance_vector_dx

  subroutine covariance_bounds(cov, x, lower, upper)
    type(covariance_model), intent(in) :: cov
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: lower(:), upper(:)
    integer :: d,j,n,pos,m
    real(dp), allocatable :: lo0(:),up0(:)
    d=size(x,2)
    if (cov%scaling) then
      allocate(lo0(d),up0(d))
      do j=1,d
        lo0(j)=1.0e-10_dp
        up0(j)=max(2.0_dp*(maxval(x(:,j))-minval(x(:,j))),1.0e-8_dp)
      end do
      n=covariance_param_count(cov,d); allocate(lower(n),upper(n)); pos=1
      do j=1,d
        m=size(cov%axis(j)%eta)
        lower(pos:pos+m-1)=1.0_dp/up0(j)
        upper(pos:pos+m-1)=1.0_dp/lo0(j)
        pos=pos+m
      end do
    else if (cov%iso) then
      allocate(lower(1),upper(1)); lower=1.0e-10_dp
      upper(1)=0.0_dp
      do j=1,d
        upper(1)=max(upper(1),2.0_dp*(maxval(x(:,j))-minval(x(:,j))))
      end do
      upper(1)=max(upper(1),1.0e-8_dp)
    else
      n=d; if(cov%kind==cov_powexp)n=2*d
      allocate(lower(n),upper(n)); lower=1.0e-10_dp
      do j=1,d
        upper(j)=max(2.0_dp*(maxval(x(:,j))-minval(x(:,j))),1.0e-8_dp)
      end do
      if(cov%kind==cov_powexp) upper(d+1:2*d)=2.0_dp
    end if
  end subroutine covariance_bounds

  pure real(dp) function scaling_fun1d(x,knots,eta) result(y)
    real(dp), intent(in) :: x,knots(:),eta(:)
    integer :: m,j
    real(dp) :: s,h,dl,dr,t
    m=size(knots)
    if(m==1) then
      y=eta(1)*(x-knots(1)); return
    end if
    if(x<knots(1)) then
      y=eta(1)*(x-knots(1)); return
    end if
    s=0.0_dp
    do j=1,m-1
      h=knots(j+1)-knots(j)
      if(x<=knots(j+1)) then
        dl=x-knots(j); dr=knots(j+1)-x; t=dl/h
        y=s+0.5_dp*t*(eta(j)*(h+dr)+eta(j+1)*dl)
        return
      end if
      s=s+0.5_dp*(eta(j)+eta(j+1))*h
    end do
    y=s+eta(m)*(x-knots(m))
  end function scaling_fun1d

  pure real(dp) function scaling_fun1d_dx(x,knots,eta) result(y)
    real(dp), intent(in) :: x,knots(:),eta(:)
    integer :: m,j
    real(dp) :: h
    m=size(knots)
    if(m==1) then
      y=eta(1); return
    end if
    if(x<knots(1)) then
      y=eta(1); return
    else if(x>knots(m)) then
      y=eta(m); return
    end if
    do j=1,m-1
      if(x<=knots(j+1)) then
        h=knots(j+1)-knots(j)
        y=eta(j)+(eta(j+1)-eta(j))*(x-knots(j))/h
        return
      end if
    end do
    y=eta(m)
  end function scaling_fun1d_dx

  subroutine scaling_grad1d(x,knots,grad)
    real(dp), intent(in) :: x(:),knots(:)
    real(dp), allocatable, intent(out) :: grad(:,:)
    integer :: i,j,m
    real(dp), allocatable :: e(:),ep(:),em(:)
    real(dp) :: h
    m=size(knots); allocate(grad(size(x),m)); allocate(e(m),ep(m),em(m)); e=1.0_dp
    do j=1,m
      h=1.0e-6_dp; ep=e; em=e; ep(j)=ep(j)+h; em(j)=em(j)-h
      do i=1,size(x)
        grad(i,j)=(scaling_fun1d(x(i),knots,ep)-scaling_fun1d(x(i),knots,em))/(2.0_dp*h)
      end do
    end do
  end subroutine scaling_grad1d

  subroutine scaling_apply(x,axis,y)
    real(dp), intent(in) :: x(:,:)
    type(scaling_axis), intent(in) :: axis(:)
    real(dp), allocatable, intent(out) :: y(:,:)
    integer :: i,j
    allocate(y(size(x,1),size(x,2)))
    do j=1,size(x,2)
      do i=1,size(x,1)
        y(i,j)=scaling_fun1d(x(i,j),axis(j)%knots,axis(j)%eta)
      end do
    end do
  end subroutine scaling_apply

end module dk_covariance
