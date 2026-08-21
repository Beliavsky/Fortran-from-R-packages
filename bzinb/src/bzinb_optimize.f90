module bzinb_optimize
  use bzinb_kinds, only : dp
  implicit none
  private
  public :: nelder_mead, golden_max
  abstract interface
    function objective_fn(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_fn
  end interface
contains
  subroutine nelder_mead(fn, x, fbest, maxiter, tol, iter, converged, step)
    procedure(objective_fn) :: fn
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: fbest
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol, step
    integer, intent(out), optional :: iter
    logical, intent(out), optional :: converged
    integer :: n, miter, it, j, ilo, ihi, inhi
    real(dp) :: ftol, st, fr, fe, fc, spread
    real(dp), allocatable :: p(:,:), f(:), c(:), xr(:), xe(:), xc(:)
    n=size(x); miter=3000; ftol=1.0e-8_dp; st=0.15_dp
    if (present(maxiter)) miter=maxiter
    if (present(tol)) ftol=tol
    if (present(step)) st=step
    allocate(p(n,n+1),f(n+1),c(n),xr(n),xe(n),xc(n))
    p(:,1)=x
    do j=2,n+1
      p(:,j)=x; p(j-1,j)=p(j-1,j)+st*(1.0_dp+abs(x(j-1)))
    end do
    do j=1,n+1; f(j)=fn(p(:,j)); end do
    do it=1,miter
      ilo=minloc(f,dim=1); ihi=maxloc(f,dim=1)
      inhi=ilo
      do j=1,n+1
        if (j==ihi) cycle
        if (inhi==ilo .or. f(j)>f(inhi)) inhi=j
      end do
      spread=maxval(abs(f-f(ilo)))
      if (spread <= ftol*(1.0_dp+abs(f(ilo)))) exit
      c=0.0_dp
      do j=1,n+1; if(j/=ihi)c=c+p(:,j); end do
      c=c/real(n,dp)
      xr=c+(c-p(:,ihi)); fr=fn(xr)
      if (fr < f(ilo)) then
        xe=c+2.0_dp*(xr-c); fe=fn(xe)
        if(fe<fr)then;p(:,ihi)=xe;f(ihi)=fe;else;p(:,ihi)=xr;f(ihi)=fr;end if
      else if (fr < f(inhi)) then
        p(:,ihi)=xr; f(ihi)=fr
      else
        if (fr < f(ihi)) then
          xc=c+0.5_dp*(xr-c)
        else
          xc=c+0.5_dp*(p(:,ihi)-c)
        end if
        fc=fn(xc)
        if(fc<min(fr,f(ihi)))then
          p(:,ihi)=xc;f(ihi)=fc
        else
          ilo=minloc(f,dim=1)
          do j=1,n+1
            if(j==ilo)cycle
            p(:,j)=p(:,ilo)+0.5_dp*(p(:,j)-p(:,ilo)); f(j)=fn(p(:,j))
          end do
        end if
      end if
    end do
    ilo=minloc(f,dim=1); x=p(:,ilo); fbest=f(ilo)
    if(present(iter))iter=min(it,miter)
    if(present(converged))converged=(it<=miter)
  end subroutine nelder_mead

  subroutine golden_max(fn1, lo, hi, xbest, fbest, tol, maxiter)
    interface
      function fn1(x) result(f)
        import dp
        real(dp), intent(in) :: x
        real(dp) :: f
      end function fn1
    end interface
    real(dp), intent(in) :: lo, hi
    real(dp), intent(out) :: xbest, fbest
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp) :: a,b,c,d,fc,fd,t,gr
    integer :: it,mi
    t=1.0e-9_dp; mi=300
    if (present(tol)) t=tol
    if (present(maxiter)) mi=maxiter
    gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp; a=lo;b=hi;c=b-gr*(b-a);d=a+gr*(b-a);fc=fn1(c);fd=fn1(d)
    do it=1,mi
      if(abs(b-a)<=t*(1.0_dp+abs(c)+abs(d)))exit
      if(fc>fd)then;b=d;d=c;fd=fc;c=b-gr*(b-a);fc=fn1(c)
      else;a=c;c=d;fc=fd;d=a+gr*(b-a);fd=fn1(d);end if
    end do
    if(fc>fd)then;xbest=c;fbest=fc;else;xbest=d;fbest=fd;end if
  end subroutine golden_max
end module bzinb_optimize
