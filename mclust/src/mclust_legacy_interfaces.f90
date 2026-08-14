! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_legacy_interfaces
  use mclust_kinds, only : dp
  implicit none
  private
  public :: ms1e, ms1v, mseii, msvii, mseei, msevi, msvei, msvvi
  public :: mseee, mseev, msvev, msvvv, msevv, msvee, mseve, msvve
  public :: hc1e, hc1v, hceii, hcvii, hceee, hcvvv

  interface
    subroutine hc1e(x,n,ic,ng,ns,nd,d)
      import dp
      integer :: n,ic(*),ng,ns,nd
      real(dp) :: x(*),d(*)
    end subroutine
    subroutine hc1v(x,n,ic,ng,ns,alpha,nd,d)
      import dp
      integer :: n,ic(*),ng,ns,nd
      real(dp) :: x(*),alpha,d(*)
    end subroutine
    subroutine hceii(x,n,p,ic,ng,ns,v,nd,d)
      import dp
      integer :: n,p,ic(*),ng,ns,nd
      real(dp) :: x(n,*),v(*),d(*)
    end subroutine
    subroutine hcvii(x,n,p,ic,ng,ns,alpha,v,nd,d)
      import dp
      integer :: n,p,ic(*),ng,ns,nd
      real(dp) :: x(n,*),alpha,v(*),d(*)
    end subroutine
    subroutine hceee(x,n,p,ic,ng,ns,io,jo,v,s,u,r)
      import dp
      integer :: n,p,ic(*),ng,ns,io(*),jo(*)
      real(dp) :: x(n,*),v(*),s(*),u(*),r(*)
    end subroutine
    subroutine hcvvv(x,n,p,ic,ng,ns,alpha,beta,v,u,s,r,nd,d)
      import dp
      integer :: n,p,ic(*),ng,ns,nd
      real(dp) :: x(n,*),alpha,beta,v(*),u(p,*),s(p,*),r(p,*),d(*)
    end subroutine
    subroutine ms1e(x,z,n,g,mu,sigsq,pro)
      import dp
      integer :: n,g
      real(dp) :: x(*),z(n,*),mu(*),sigsq,pro(*)
    end subroutine
    subroutine ms1v(x,z,n,g,mu,sigsq,pro)
      import dp
      integer :: n,g
      real(dp) :: x(*),z(n,*),mu(*),sigsq(*),pro(*)
    end subroutine
    subroutine mseii(x,z,n,p,g,mu,sigsq,pro)
      import dp
      integer :: n,p,g
      real(dp) :: x(n,*),z(n,*),mu(p,*),sigsq,pro(*)
    end subroutine
    subroutine msvii(x,z,n,p,g,mu,sigsq,pro)
      import dp
      integer :: n,p,g
      real(dp) :: x(n,*),z(n,*),mu(p,*),sigsq(*),pro(*)
    end subroutine
    subroutine mseei(x,z,n,p,g,mu,scale,shape,pro)
      import dp
      integer :: n,p,g
      real(dp) :: x(n,*),z(n,*),mu(p,*),scale,shape(*),pro(*)
    end subroutine
    subroutine msevi(x,z,n,p,g,mu,scale,shape,pro)
      import dp
      integer :: n,p,g
      real(dp) :: x(n,*),z(n,*),mu(p,*),scale,shape(p,*),pro(*)
    end subroutine
    subroutine msvei(x,z,n,p,g,maxi,tol,mu,scale,shape,pro,scl,shp,w)
      import dp
      integer :: n,p,g,maxi
      real(dp) :: x(n,*),z(n,*),tol,mu(p,*),scale(*),shape(*),pro(*),scl(*),shp(*),w(p,*)
    end subroutine
    subroutine msvvi(x,z,n,p,g,mu,scale,shape,pro)
      import dp
      integer :: n,p,g
      real(dp) :: x(n,*),z(n,*),mu(p,*),scale(*),shape(p,*),pro(*)
    end subroutine
    subroutine mseee(x,z,n,p,g,w,mu,u,pro)
      import dp
      integer :: n,p,g
      real(dp) :: x(n,*),z(n,*),w(*),mu(p,*),u(p,*),pro(*)
    end subroutine
    subroutine mseev(x,z,n,p,g,w,lwork,mu,scale,shape,o,pro)
      import dp
      integer :: n,p,g,lwork
      real(dp) :: x(n,*),z(n,*),w(*),mu(p,*),scale,shape(*),o(p,p,*),pro(*)
    end subroutine
    subroutine msvev(x,z,n,p,g,w,lwork,maxi,tol,mu,scale,shape,o,pro)
      import dp
      integer :: n,p,g,lwork,maxi
      real(dp) :: x(n,*),z(n,*),w(*),tol,mu(p,*),scale(*),shape(*),o(p,p,*),pro(*)
    end subroutine
    subroutine msvvv(x,z,n,p,g,w,mu,u,pro,s)
      import dp
      integer :: n,p,g
      real(dp) :: x(n,*),z(n,*),w(*),mu(p,*),u(p,p,*),pro(*),s(p,*)
    end subroutine
    subroutine msevv(x,z,n,p,g,mu,o,u,scale,shape,pro,lwork,info,eps)
      import dp
      integer :: n,p,g,lwork,info
      real(dp) :: x(n,p),z(n,g),mu(p,g),o(p,p,*),u(p,p,*),scale(g),shape(p,g),pro(g),eps
    end subroutine
    subroutine msvee(x,z,n,p,g,mu,u,c,scale,pro,lwork,info,itmax,tol,niterin,errin)
      import dp
      integer :: n,p,g,lwork,info,itmax,niterin
      real(dp) :: x(n,p),z(n,g),mu(p,g),u(p,p,g),c(p,p),scale(g),pro(g),tol,errin
    end subroutine
    subroutine mseve(x,z,n,p,g,mu,u,o,scale,shape,pro,lwork,info,itmax,tol,niterin,errin,eps)
      import dp
      integer :: n,p,g,lwork,info,itmax,niterin
      real(dp) :: x(n,p),z(n,g),mu(p,g),u(p,p,g),o(p,p),scale,shape(p,g),pro(g),tol,errin,eps
    end subroutine
    subroutine msvve(x,z,n,p,g,mu,u,o,scale,shape,pro,lwork,info,itmax,tol,niterin,errin,eps)
      import dp
      integer :: n,p,g,lwork,info,itmax,niterin
      real(dp) :: x(n,p),z(n,g),mu(p,g),u(p,p,g),o(p,p),scale(g),shape(p,g),pro(g),tol,errin,eps
    end subroutine
  end interface
end module mclust_legacy_interfaces
