! SPDX-License-Identifier: GPL-2.0-or-later
subroutine desolve_zvode_entry(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt, &
    zwork,lzw,rwork,lrw,iwork,liw,jac,mf,rpar,ipar)
  implicit none
  external :: f,jac,zvode
  integer :: neq,itol,itask,istate,iopt,lzw,lrw,iwork(*),liw,mf,ipar(*)
  complex(kind=kind(1.0d0)) :: y(*),zwork(*)
  double precision :: t,tout,rtol(*),atol(*),rwork(*),rpar(*)
  call zvode(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,zwork,lzw,rwork, &
      lrw,iwork,liw,jac,mf,rpar,ipar)
end subroutine desolve_zvode_entry
