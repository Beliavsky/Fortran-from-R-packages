! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! This translation is distributed under the same license choice; see
! LICENSE-GPL-2 and LICENSE-GPL-3 in the project root.
module dk_optimizer
  use dk_kinds, only : dp
  implicit none
  private

  abstract interface
    function objective_iface(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_iface
    subroutine gradient_iface(x,g)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
    end subroutine gradient_iface
  end interface

  public :: bounded_bfgs

contains

  subroutine bounded_bfgs(fun, x0, lower, upper, xbest, fbest, info, gradient, max_iter, tol)
    procedure(objective_iface) :: fun
    real(dp), intent(in) :: x0(:), lower(:), upper(:)
    real(dp), allocatable, intent(out) :: xbest(:)
    real(dp), intent(out) :: fbest
    integer, intent(out) :: info
    procedure(gradient_iface), optional :: gradient
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol

    integer :: n, iter, mit, ls, i
    real(dp) :: epsg, f, fn, alpha, slope, ys, rho
    real(dp), allocatable :: z(:), zn(:), x(:), xn(:), g(:), gn(:), gz(:), gzn(:)
    real(dp), allocatable :: h(:,:), p(:), s(:), y(:), v(:,:), ident(:,:)

    n=size(x0); mit=300; if(present(max_iter))mit=max_iter
    epsg=1.0e-7_dp; if(present(tol))epsg=tol
    allocate(z(n),zn(n),x(n),xn(n),g(n),gn(n),gz(n),gzn(n),h(n,n),p(n),s(n),y(n),v(n,n),ident(n,n))
    ident=0.0_dp
    do i=1,n
      ident(i,i)=1.0_dp
    end do
    call x_to_z(x0,lower,upper,z)
    call z_to_x(z,lower,upper,x)
    f=fun(x)
    call grad_x(fun,x,lower,upper,g,gradient)
    call transform_grad(z,lower,upper,g,gz)
    h=ident
    info=1

    do iter=1,mit
      if(maxval(abs(gz)) <= epsg) then
        info=0; exit
      end if
      p=-matmul(h,gz)
      slope=dot_product(gz,p)
      if(slope >= -epsilon(1.0_dp)) then
        p=-gz; h=ident; slope=-dot_product(gz,gz)
      end if
      alpha=1.0_dp
      do ls=1,35
        zn=z+alpha*p
        call z_to_x(zn,lower,upper,xn)
        fn=fun(xn)
        if(fn <= f+1.0e-4_dp*alpha*slope) exit
        alpha=0.5_dp*alpha
      end do
      if(ls>35) exit
      call grad_x(fun,xn,lower,upper,gn,gradient)
      call transform_grad(zn,lower,upper,gn,gzn)
      s=zn-z; y=gzn-gz; ys=dot_product(y,s)
      if(ys > 1.0e-12_dp*sqrt(max(dot_product(y,y)*dot_product(s,s),tiny(1.0_dp)))) then
        rho=1.0_dp/ys
        v=ident-rho*outer(s,y)
        h=matmul(v,matmul(h,transpose(v)))+rho*outer(s,s)
      else
        h=ident
      end if
      z=zn; x=xn; f=fn; g=gn; gz=gzn
    end do
    xbest=x; fbest=f
  contains
    function outer(a,b) result(c)
      real(dp), intent(in) :: a(:),b(:)
      real(dp) :: c(size(a),size(b))
      integer :: ii,jj
      do jj=1,size(b)
        do ii=1,size(a)
          c(ii,jj)=a(ii)*b(jj)
        end do
      end do
    end function outer
  end subroutine bounded_bfgs

  subroutine z_to_x(z,lo,hi,x)
    real(dp), intent(in) :: z(:),lo(:),hi(:)
    real(dp), intent(out) :: x(:)
    integer :: i
    real(dp) :: s
    do i=1,size(z)
      if(z(i)>=0.0_dp) then
        s=1.0_dp/(1.0_dp+exp(-min(z(i),700.0_dp)))
      else
        s=exp(max(z(i),-700.0_dp))/(1.0_dp+exp(max(z(i),-700.0_dp)))
      end if
      x(i)=lo(i)+(hi(i)-lo(i))*s
    end do
  end subroutine z_to_x

  subroutine x_to_z(x,lo,hi,z)
    real(dp), intent(in) :: x(:),lo(:),hi(:)
    real(dp), intent(out) :: z(:)
    integer :: i
    real(dp) :: q
    do i=1,size(x)
      q=(x(i)-lo(i))/(hi(i)-lo(i))
      q=min(max(q,1.0e-15_dp),1.0_dp-1.0e-15_dp)
      z(i)=log(q/(1.0_dp-q))
    end do
  end subroutine x_to_z

  subroutine transform_grad(z,lo,hi,g,gz)
    real(dp), intent(in) :: z(:),lo(:),hi(:),g(:)
    real(dp), intent(out) :: gz(:)
    integer :: i
    real(dp) :: s
    do i=1,size(z)
      if(z(i)>=0.0_dp) then
        s=1.0_dp/(1.0_dp+exp(-min(z(i),700.0_dp)))
      else
        s=exp(max(z(i),-700.0_dp))/(1.0_dp+exp(max(z(i),-700.0_dp)))
      end if
      gz(i)=g(i)*(hi(i)-lo(i))*s*(1.0_dp-s)
    end do
  end subroutine transform_grad

  subroutine grad_x(fun,x,lo,hi,g,gradient)
    procedure(objective_iface) :: fun
    real(dp), intent(in) :: x(:),lo(:),hi(:)
    real(dp), intent(out) :: g(:)
    procedure(gradient_iface), optional :: gradient
    real(dp), allocatable :: xp(:),xm(:)
    real(dp) :: h,fp,fm
    integer :: i
    if(present(gradient)) then
      call gradient(x,g)
      return
    end if
    xp=x; xm=x
    do i=1,size(x)
      h=max(1.0e-7_dp,1.0e-5_dp*max(1.0_dp,abs(x(i))))
      xp=x; xm=x
      xp(i)=min(hi(i),x(i)+h); xm(i)=max(lo(i),x(i)-h)
      fp=fun(xp); fm=fun(xm)
      g(i)=(fp-fm)/(xp(i)-xm(i))
    end do
  end subroutine grad_x

end module dk_optimizer
