! SPDX-License-Identifier: GPL-2.0-only
module calibrar_stopping
  use calibrar_kinds, only : dp
  implicit none
  private
  public :: smooth_stop2, smooth_stop3, smooth_stop4, n_stop
contains
  function smooth_stop2(x,reltol,nwin) result(stop)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::reltol
    integer,intent(in),optional::nwin
    logical::stop
    real(dp)::rt,x0,x1,dx,tol
    integer::n,w
    rt=sqrt(epsilon(1.0_dp));if(present(reltol))rt=reltol
    w=10;if(present(nwin))w=nwin;n=size(x);stop=.false.;if(n<3*w)return
    x0=sum(x(n-3*w+1:n-w))/real(2*w,dp)
    x1=sum(x(n-2*w+1:n))/real(2*w,dp)
    dx=x0-x1;tol=rt*(abs(x0)+rt);stop=abs(dx)<tol
  end function smooth_stop2

  function smooth_stop3(x,reltol,nwin) result(stop)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::reltol
    integer,intent(in),optional::nwin
    logical::stop
    real(dp)::rt,x0,x1,dx,tol
    integer::n,w
    rt=sqrt(epsilon(1.0_dp));if(present(reltol))rt=reltol
    w=10;if(present(nwin))w=nwin;n=size(x);stop=.false.;if(n<2*w)return
    x0=sum(x(n-2*w+1:n-w))/real(w,dp)
    x1=sum(x(n-w+1:n))/real(w,dp)
    dx=x0-x1;tol=rt*(abs(x0)+rt);stop=abs(dx)<tol
  end function smooth_stop3

  function smooth_stop4(x,reltol,nwin) result(stop)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::reltol
    integer,intent(in),optional::nwin
    logical::stop
    real(dp)::rt,x0,x1,dx,tol
    integer::n,w
    rt=sqrt(epsilon(1.0_dp));if(present(reltol))rt=reltol
    w=10;if(present(nwin))w=nwin;n=size(x);stop=.false.;if(n<2*w)return
    rt=real(w,dp)*rt
    x0=maxval(x(n-w+1:n));x1=minval(x(n-w+1:n));dx=x0-x1
    tol=rt*(abs(x0)+rt);stop=dx<tol
  end function smooth_stop4

  function n_stop(x,nwin) result(stop)
    logical,intent(in)::x(:)
    integer,intent(in)::nwin
    logical::stop
    integer::n,i,count
    n=size(x);stop=.false.;if(n<=2*nwin)return
    count=0
    do i=n-nwin+1,n
      if(x(i))count=count+1
    end do
    stop=real(count,dp)/real(nwin,dp)>(1.0_dp-0.1_dp/real(nwin,dp))
  end function n_stop
end module calibrar_stopping
