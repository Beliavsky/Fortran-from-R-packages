module arfima_optimizer
  use arfima_kinds, only : dp
  use arfima_status, only : arfima_ok, arfima_invalid_input, arfima_no_convergence
  implicit none
  private
  public :: objective_function, nelder_mead, numerical_hessian

  abstract interface
    function objective_function(x,context) result(value)
      import :: dp
      real(dp),intent(in)::x(:)
      class(*),intent(inout)::context
      real(dp)::value
    end function objective_function
  end interface

contains

  subroutine nelder_mead(fn,context,x,fval,iterations,evaluations,info,max_iterations,tolerance,initial_step)
    procedure(objective_function)::fn
    class(*),intent(inout)::context
    real(dp),intent(inout)::x(:)
    real(dp),intent(out)::fval
    integer,intent(out)::iterations,evaluations,info
    integer,intent(in),optional::max_iterations
    real(dp),intent(in),optional::tolerance,initial_step
    real(dp),allocatable::simplex(:,:),values(:),centroid(:),xr(:),xe(:),xc(:)
    real(dp)::alpha,gamma,rho,sigma,tol,step,fr,fe,fc,spread,size_simplex
    integer::n,maxit,j,best,worst,second_worst

    n=size(x); maxit=2000; if(present(max_iterations)) maxit=max_iterations
    tol=1.0e-7_dp; if(present(tolerance)) tol=tolerance
    step=0.15_dp; if(present(initial_step)) step=initial_step
    alpha=1.0_dp; gamma=2.0_dp; rho=0.5_dp; sigma=0.5_dp
    if(n==0) then
      fval=fn(x,context); iterations=0; evaluations=1; info=arfima_ok; return
    end if
    allocate(simplex(n,n+1),values(n+1),centroid(n),xr(n),xe(n),xc(n))
    simplex(:,1)=x
    do j=1,n
      simplex(:,j+1)=x
      simplex(j,j+1)=simplex(j,j+1)+step*max(1.0_dp,abs(x(j)))
    end do
    do j=1,n+1; values(j)=fn(simplex(:,j),context); end do
    evaluations=n+1; info=arfima_no_convergence
    do iterations=1,maxit
      call order_simplex(values,best,worst,second_worst)
      spread=maxval(abs(values-values(best)))
      size_simplex=0.0_dp
      do j=1,n+1; size_simplex=max(size_simplex,maxval(abs(simplex(:,j)-simplex(:,best)))); end do
      if(spread<=tol*(1.0_dp+abs(values(best))) .and. size_simplex<=sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,best))))) then
        info=arfima_ok; exit
      end if
      centroid=0.0_dp
      do j=1,n+1; if(j/=worst) centroid=centroid+simplex(:,j); end do
      centroid=centroid/real(n,dp)
      xr=centroid+alpha*(centroid-simplex(:,worst)); fr=fn(xr,context); evaluations=evaluations+1
      if(fr<values(best)) then
        xe=centroid+gamma*(xr-centroid); fe=fn(xe,context); evaluations=evaluations+1
        if(fe<fr) then; simplex(:,worst)=xe; values(worst)=fe
        else; simplex(:,worst)=xr; values(worst)=fr; end if
      else if(fr<values(second_worst)) then
        simplex(:,worst)=xr; values(worst)=fr
      else
        if(fr<values(worst)) then
          xc=centroid+rho*(xr-centroid)
        else
          xc=centroid-rho*(centroid-simplex(:,worst))
        end if
        fc=fn(xc,context); evaluations=evaluations+1
        if(fc<min(fr,values(worst))) then
          simplex(:,worst)=xc; values(worst)=fc
        else
          do j=1,n+1
            if(j/=best) then
              simplex(:,j)=simplex(:,best)+sigma*(simplex(:,j)-simplex(:,best))
              values(j)=fn(simplex(:,j),context)
            end if
          end do
          evaluations=evaluations+n
        end if
      end if
    end do
    call order_simplex(values,best,worst,second_worst)
    x=simplex(:,best); fval=values(best)
  contains
    subroutine order_simplex(v,ib,iw,isw)
      real(dp),intent(in)::v(:)
      integer,intent(out)::ib,iw,isw
      integer::k
      ib=minloc(v,dim=1); iw=maxloc(v,dim=1); isw=ib
      do k=1,size(v)
        if(k/=iw) then
          if(isw==iw .or. v(k)>v(isw)) isw=k
        end if
      end do
    end subroutine order_simplex
  end subroutine nelder_mead

  subroutine numerical_hessian(fn,context,x,hessian,info,relative_step)
    procedure(objective_function)::fn
    class(*),intent(inout)::context
    real(dp),intent(in)::x(:)
    real(dp),allocatable,intent(out)::hessian(:,:)
    integer,intent(out)::info
    real(dp),intent(in),optional::relative_step
    real(dp),allocatable::xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:),h(:)
    real(dp)::f0,fp,fm,fpp,fpm,fmp,fmm,rel
    integer::n,i,j
    n=size(x); rel=epsilon(1.0_dp)**0.25_dp; if(present(relative_step)) rel=relative_step
    allocate(hessian(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n),h(n))
    h=rel*max(1.0_dp,abs(x)); hessian=0.0_dp; f0=fn(x,context)
    if(.not.(f0<huge(1.0_dp))) then; info=arfima_invalid_input; return; end if
    do i=1,n
      xp=x; xm=x; xp(i)=xp(i)+h(i); xm(i)=xm(i)-h(i)
      fp=fn(xp,context); fm=fn(xm,context)
      hessian(i,i)=(fp-2.0_dp*f0+fm)/(h(i)*h(i))
      do j=i+1,n
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=xpp(i)+h(i); xpp(j)=xpp(j)+h(j)
        xpm(i)=xpm(i)+h(i); xpm(j)=xpm(j)-h(j)
        xmp(i)=xmp(i)-h(i); xmp(j)=xmp(j)+h(j)
        xmm(i)=xmm(i)-h(i); xmm(j)=xmm(j)-h(j)
        fpp=fn(xpp,context); fpm=fn(xpm,context); fmp=fn(xmp,context); fmm=fn(xmm,context)
        hessian(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*h(i)*h(j)); hessian(j,i)=hessian(i,j)
      end do
    end do
    info=arfima_ok
  end subroutine numerical_hessian
end module arfima_optimizer
