! SPDX-License-Identifier: BSD-3-Clause
module rfortran_arpack
  implicit none
  private
  public :: dsaupd,dseupd,dnaupd,dneupd

  interface
    subroutine dsaupd(ido,bm,n,wh,nev,tol,res,ncv,v,ldv,ip,ipt,wd,wl,lwl,info)
      integer :: ido,n,nev,ncv,ldv,lwl,info,ip(*),ipt(*)
      character(len=1) :: bm
      character(len=2) :: wh
      double precision :: tol,res(*),v(ldv,*),wd(*),wl(*)
    end subroutine dsaupd

    subroutine dseupd(rv,hm,sel,d,z,ldz,sig,bm,n,wh,nev,tol,res,ncv,v,ldv,ip,ipt,wd,wl,lwl,info)
      logical :: rv,sel(*)
      character(len=1) :: hm,bm
      character(len=2) :: wh
      integer :: ldz,n,nev,ncv,ldv,lwl,info,ip(*),ipt(*)
      double precision :: d(*),z(ldz,*),sig,tol,res(*),v(ldv,*),wd(*),wl(*)
    end subroutine dseupd

    subroutine dnaupd(ido,bm,n,wh,nev,tol,res,ncv,v,ldv,ip,ipt,wd,wl,lwl,info)
      integer :: ido,n,nev,ncv,ldv,lwl,info,ip(*),ipt(*)
      character(len=1) :: bm
      character(len=2) :: wh
      double precision :: tol,res(*),v(ldv,*),wd(*),wl(*)
    end subroutine dnaupd

    subroutine dneupd(rv,hm,sel,dr,di,z,ldz,sr,si,we,bm,n,wh,nev,tol,res,ncv,v,ldv,ip,ipt,wd,wl,lwl,info)
      logical :: rv,sel(*)
      character(len=1) :: hm,bm
      character(len=2) :: wh
      integer :: ldz,n,nev,ncv,ldv,lwl,info,ip(*),ipt(*)
      double precision :: dr(*),di(*),z(ldz,*),sr,si,we(*),tol,res(*),v(ldv,*),wd(*),wl(*)
    end subroutine dneupd
  end interface
end module rfortran_arpack
