! SPDX-License-Identifier: GPL-2.0-or-later
module desolve_callback_bridge
  use desolve_kinds, only : dp
  use desolve_types, only : ode_rhs, complex_ode_rhs, ode_root, dae_residual
  implicit none
  private
  procedure(ode_rhs), pointer, save, public :: current_rhs => null()
  procedure(complex_ode_rhs), pointer, save, public :: current_complex_rhs => null()
  procedure(ode_root), pointer, save, public :: current_root => null()
  procedure(dae_residual), pointer, save, public :: current_residual => null()
  integer, save, public :: current_n = 0
  integer, save, public :: current_ng = 0
  public :: odepack_rhs_bridge, odepack_jac_dummy, odepack_root_bridge, lsodes_jac_dummy
  public :: vode_rhs_bridge, vode_jac_dummy
  public :: zvode_rhs_bridge, zvode_jac_dummy
  public :: daspk_res_bridge, daspk_jac_dummy, daspk_psol_dummy
  public :: radau_rhs_bridge, radau_jac_dummy, radau_mass_identity, radau_solout_dummy
contains

  subroutine odepack_rhs_bridge(neq, t, y, ydot, rpar, ipar)
    integer, intent(in) :: neq
    real(dp), intent(in) :: t, y(*)
    real(dp), intent(out) :: ydot(*)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    real(dp) :: yin(neq), dy(neq)
    integer :: i
    if (.not. associated(current_rhs)) error stop 'deSolve: ODE RHS callback not set'
    do i = 1, neq
      yin(i) = y(i)
    end do
    call current_rhs(t, yin, dy)
    do i = 1, neq
      ydot(i) = dy(i)
    end do
    if (size_hint(rpar, ipar) < 0) error stop 'unreachable'
  end subroutine odepack_rhs_bridge

  subroutine vode_rhs_bridge(neq, t, y, ydot, rpar, ipar)
    integer, intent(in) :: neq
    real(dp), intent(in) :: t, y(*)
    real(dp), intent(out) :: ydot(*)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    call odepack_rhs_bridge(neq, t, y, ydot, rpar, ipar)
  end subroutine vode_rhs_bridge

  subroutine odepack_jac_dummy(neq, t, y, ml, mu, pd, nrowpd, rpar, ipar)
    integer, intent(in) :: neq, ml, mu, nrowpd
    real(dp), intent(in) :: t, y(*)
    real(dp), intent(out) :: pd(nrowpd,*)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    integer :: i, j
    do j = 1, neq
      do i = 1, nrowpd
        pd(i,j) = 0.0_dp
      end do
    end do
    if (t + y(1) + real(ml+mu,dp) + rpar(1) + real(ipar(1),dp) < -huge(1.0_dp)) stop
  end subroutine odepack_jac_dummy

  subroutine vode_jac_dummy(neq, t, y, ml, mu, pd, nrowpd, rpar, ipar)
    integer, intent(in) :: neq, ml, mu, nrowpd
    real(dp), intent(in) :: t, y(*)
    real(dp), intent(out) :: pd(nrowpd,*)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    call odepack_jac_dummy(neq,t,y,ml,mu,pd,nrowpd,rpar,ipar)
  end subroutine vode_jac_dummy

  subroutine lsodes_jac_dummy(neq, t, y, jcol, ian, jan, pdj, rpar, ipar)
    integer, intent(in) :: neq, jcol
    real(dp), intent(in) :: t, y(*)
    integer, intent(in) :: ian(*), jan(*)
    real(dp), intent(out) :: pdj(*)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    integer :: i
    do i=1,neq; pdj(i)=0.0_dp; end do
    if (t+y(1)+real(jcol+ian(1)+jan(1)+ipar(1),dp)+rpar(1) < -huge(1.0_dp)) stop
  end subroutine lsodes_jac_dummy

  subroutine odepack_root_bridge(neq, t, y, ng, gout, rpar, ipar)
    integer, intent(in) :: neq, ng
    real(dp), intent(in) :: t, y(*)
    real(dp), intent(out) :: gout(*)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    real(dp) :: yin(neq), g(ng)
    integer :: i
    if (.not. associated(current_root)) error stop 'deSolve: root callback not set'
    do i=1,neq; yin(i)=y(i); end do
    call current_root(t, yin, g)
    do i=1,ng; gout(i)=g(i); end do
    if (rpar(1) + real(ipar(1),dp) < -huge(1.0_dp)) stop
  end subroutine odepack_root_bridge

  subroutine zvode_rhs_bridge(neq, t, y, ydot, rpar, ipar)
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    complex(dp), intent(in) :: y(*)
    real(dp), intent(inout) :: rpar(*)
    complex(dp), intent(out) :: ydot(*)
    integer, intent(inout) :: ipar(*)
    complex(dp) :: yin(neq), dy(neq)
    integer :: i
    if (.not. associated(current_complex_rhs)) error stop 'deSolve: complex RHS callback not set'
    do i=1,neq; yin(i)=y(i); end do
    call current_complex_rhs(t, yin, dy)
    do i=1,neq; ydot(i)=dy(i); end do
    if (rpar(1) + real(ipar(1),dp) < -huge(1.0_dp)) stop
  end subroutine zvode_rhs_bridge

  subroutine zvode_jac_dummy(neq, t, y, ml, mu, pd, nrowpd, rpar, ipar)
    integer, intent(in) :: neq, ml, mu, nrowpd
    real(dp), intent(in) :: t
    complex(dp), intent(in) :: y(*)
    real(dp), intent(inout) :: rpar(*)
    complex(dp), intent(out) :: pd(nrowpd,*)
    integer, intent(inout) :: ipar(*)
    integer :: i,j
    do j=1,neq
      do i=1,nrowpd
        pd(i,j)=(0.0_dp,0.0_dp)
      end do
    end do
    if (t + real(y(1),dp) + rpar(1) + real(ml+mu+ipar(1),dp) < -huge(1.0_dp)) stop
  end subroutine zvode_jac_dummy

  subroutine daspk_res_bridge(t, y, yprime, cj, delta, ires, rpar, ipar)
    real(dp), intent(in) :: t, y(*), yprime(*), cj
    real(dp), intent(out) :: delta(*)
    integer, intent(inout) :: ires, ipar(*)
    real(dp), intent(inout) :: rpar(*)
    real(dp) :: yy(current_n), yp(current_n), dd(current_n)
    integer :: i
    if (.not. associated(current_residual)) error stop 'deSolve: DAE residual callback not set'
    do i=1,current_n; yy(i)=y(i); yp(i)=yprime(i); end do
    call current_residual(t,yy,yp,dd)
    do i=1,current_n; delta(i)=dd(i); end do
    ires = 0
    if (cj + rpar(1) + real(ipar(1),dp) < -huge(1.0_dp)) stop
  end subroutine daspk_res_bridge

  subroutine daspk_jac_dummy(t, y, yprime, pd, cj, rpar, ipar)
    real(dp), intent(in) :: t, y(*), yprime(*), cj
    real(dp), intent(out) :: pd(*)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    integer :: i
    do i=1,current_n*current_n; pd(i)=0.0_dp; end do
    if (t+y(1)+yprime(1)+cj+rpar(1)+real(ipar(1),dp) < -huge(1.0_dp)) stop
  end subroutine daspk_jac_dummy

  subroutine daspk_psol_dummy(neq, t, y, yprime, savr, wk, cj, wght, wp, iwp, &
      b, eplin, ier, rpar, ipar)
    integer, intent(in) :: neq
    real(dp), intent(in) :: t,y(*),yprime(*),savr(*),wk(*),cj,wght(*),wp(*),eplin
    integer, intent(in) :: iwp(*)
    real(dp), intent(inout) :: b(*), rpar(*)
    integer, intent(out) :: ier
    integer, intent(inout) :: ipar(*)
    ier=0
    if (neq + iwp(1) + ipar(1) < -huge(1)) stop
    if (t+y(1)+yprime(1)+savr(1)+wk(1)+cj+wght(1)+wp(1)+b(1)+eplin+rpar(1) < -huge(1.0_dp)) stop
  end subroutine daspk_psol_dummy

  subroutine radau_rhs_bridge(n, x, y, f, rpar, ipar)
    integer, intent(in) :: n
    real(dp), intent(in) :: x,y(n)
    real(dp), intent(out) :: f(n)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    if (.not. associated(current_rhs)) error stop 'deSolve: RADAU RHS callback not set'
    call current_rhs(x,y,f)
    if (rpar(1)+real(ipar(1),dp) < -huge(1.0_dp)) stop
  end subroutine radau_rhs_bridge

  subroutine radau_jac_dummy(n, x, y, dfy, ldfy, rpar, ipar)
    integer, intent(in) :: n, ldfy
    real(dp), intent(in) :: x,y(n)
    real(dp), intent(out) :: dfy(ldfy,n)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    dfy=0.0_dp
    if (x+y(1)+rpar(1)+real(ipar(1),dp) < -huge(1.0_dp)) stop
  end subroutine radau_jac_dummy

  subroutine radau_mass_identity(n, am, lmas, rpar, ipar)
    integer, intent(in) :: n, lmas
    real(dp), intent(out) :: am(lmas,n)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    integer :: i
    am=0.0_dp
    do i=1,min(n,lmas); am(i,i)=1.0_dp; end do
    if (rpar(1)+real(ipar(1),dp) < -huge(1.0_dp)) stop
  end subroutine radau_mass_identity

  subroutine radau_solout_dummy(nr, xold, x, y, cont, lrc, n, rpar, ipar, irtrn)
    integer, intent(in) :: nr,lrc,n
    real(dp), intent(in) :: xold,x,y(n),cont(lrc)
    real(dp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*), irtrn
    irtrn=1
    if (nr+lrc+n+ipar(1) < -huge(1)) stop
    if (xold+x+y(1)+cont(1)+rpar(1) < -huge(1.0_dp)) stop
  end subroutine radau_solout_dummy

  integer function size_hint(rpar, ipar) result(v)
    real(dp), intent(in) :: rpar(*)
    integer, intent(in) :: ipar(*)
    v=0
    if (rpar(1)+real(ipar(1),dp) < -huge(1.0_dp)) v=-1
  end function size_hint

end module desolve_callback_bridge
