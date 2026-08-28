! SPDX-License-Identifier: GPL-2.0-or-later
module nleqslv_fortran
   use, intrinsic :: iso_fortran_env, only : real64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use nleqslv_legacy_kernels, only : liqsiz, nwnleq
   implicit none
   private

   integer, parameter, public :: dp = real64
   integer, parameter, public :: NLEQ_NEWTON = 0, NLEQ_BROYDEN = 1
   integer, parameter, public :: NLEQ_NONE = 0, NLEQ_CLINE = 1, NLEQ_QLINE = 2
   integer, parameter, public :: NLEQ_GLINE = 3, NLEQ_DBLDOG = 4, NLEQ_PWLDOG = 5, NLEQ_HOOK = 6
   integer, parameter, public :: NLEQ_SCALE_FIXED = 0, NLEQ_SCALE_AUTO = 1

   type, public :: nleq_options
      integer :: method = NLEQ_BROYDEN
      integer :: global = NLEQ_DBLDOG
      integer :: xscalm = NLEQ_SCALE_FIXED
      integer :: maxit = 150
      real(dp) :: ftol = 1.0e-8_dp
      real(dp) :: xtol = 1.0e-8_dp
      real(dp) :: btol = 1.0e-3_dp
      real(dp) :: cndtol = 1.0e-12_dp
      real(dp) :: stepmax = -1.0_dp
      real(dp) :: delta = -2.0_dp
      real(dp) :: sigma = 0.5_dp
      logical :: allow_singular = .false.
      logical :: check_jacobian = .false.
      logical :: return_jacobian = .false.
      logical :: trace = .false.
      integer :: dsub = -1
      integer :: dsuper = -1
      real(dp), allocatable :: scalex(:)
   end type nleq_options

   type, public :: nleq_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: fvec(:)
      real(dp), allocatable :: scalex(:)
      real(dp), allocatable :: jac(:,:)
      integer :: termcd = 0
      integer :: nfcnt = 0
      integer :: njcnt = 0
      integer :: iter = 0
      real(dp) :: jac_rcond = 0.0_dp
      character(len=160) :: message = ''
   end type nleq_result

   type, public :: search_zeros_result
      real(dp), allocatable :: x(:,:)
      real(dp), allocatable :: xfnorm(:)
      integer, allocatable :: converged_rows(:)
      integer, allocatable :: xtol_rows(:)
      integer, allocatable :: failed_rows(:)
   end type search_zeros_result

   type, public :: nleq_test_result
      integer, allocatable :: method(:)
      integer, allocatable :: global(:)
      integer, allocatable :: termcd(:)
      integer, allocatable :: nfcnt(:)
      integer, allocatable :: njcnt(:)
      integer, allocatable :: iter(:)
      real(dp), allocatable :: fnorm(:)
      real(dp), allocatable :: cpu_seconds(:)
   end type nleq_test_result

   abstract interface
      subroutine nleq_function(x, f)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: f(:)
      end subroutine nleq_function
      subroutine nleq_jacobian(x, jac)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: jac(:,:)
      end subroutine nleq_jacobian
   end interface

   procedure(nleq_function), pointer, save :: active_fn => null()
   procedure(nleq_jacobian), pointer, save :: active_jac => null()
   logical, save :: solver_active = .false.

   public :: solve_nleqslv, search_zeros, test_nleqslv, termination_message
   public :: nleq_function, nleq_jacobian


contains

   subroutine solve_nleqslv(x0, fn, result, options, jac)
      real(dp), intent(in) :: x0(:)
      procedure(nleq_function) :: fn
      type(nleq_result), intent(out) :: result
      type(nleq_options), intent(in), optional :: options
      procedure(nleq_jacobian), optional :: jac

      type(nleq_options) :: opt
      integer :: n, lrwork, qrwsiz, ldr, lrjac
      integer :: maxit, method, global, xscalm, njcnt, nfcnt, iter, termcd
      integer :: jacflg(4), outopt(3)
      real(dp) :: xtol, ftol, btol, cndtol, stepmx, delta, sigma
      real(dp), allocatable :: x(:), xp(:), fp(:), gp(:), scalex(:)
      real(dp), allocatable :: rwork(:), rcdwrk(:), qrwork(:), rjac(:,:)
      integer, allocatable :: icdwrk(:)
      real(dp) :: fcheck(size(x0))

      if (solver_active) error stop 'solve_nleqslv: recursive/concurrent call is not supported'
      n = size(x0)
      if (n < 1) error stop 'solve_nleqslv: x0 must be nonempty'
      if (any(.not. ieee_is_finite(x0))) error stop 'solve_nleqslv: x0 contains non-finite values'

      opt = nleq_options()
      if (present(options)) opt = options
      if (opt%global == NLEQ_NONE .and. opt%maxit == 150) opt%maxit = 20

      call fn(x0, fcheck)
      if (any(.not. ieee_is_finite(fcheck))) error stop 'solve_nleqslv: initial function value is non-finite'

      method = opt%method; global = opt%global; xscalm = opt%xscalm; maxit = opt%maxit
      xtol = opt%xtol; ftol = opt%ftol; btol = opt%btol; cndtol = opt%cndtol
      stepmx = opt%stepmax; delta = opt%delta; sigma = opt%sigma
      if (any([.not.ieee_is_finite(xtol), .not.ieee_is_finite(ftol), .not.ieee_is_finite(btol), &
               .not.ieee_is_finite(cndtol), .not.ieee_is_finite(stepmx), .not.ieee_is_finite(delta), &
               .not.ieee_is_finite(sigma)])) error stop 'solve_nleqslv: non-finite control value'

      jacflg = 0
      if (present(jac)) jacflg(1) = 1
      if (opt%dsub >= 0 .or. opt%dsuper >= 0) then
         if (opt%dsub < 0 .or. opt%dsuper < 0) error stop 'solve_nleqslv: dsub and dsuper must both be specified'
         if (opt%dsub > n-2 .or. opt%dsuper > n-2) error stop 'solve_nleqslv: invalid band widths'
         jacflg(1) = jacflg(1) + 2
         jacflg(2) = opt%dsub
         jacflg(3) = opt%dsuper
      else
         jacflg(2:3) = -1
      end if
      jacflg(4) = merge(1, 0, opt%allow_singular)
      if (method == NLEQ_BROYDEN .and. opt%dsub == 0 .and. opt%dsuper == 0) &
         error stop 'solve_nleqslv: Broyden is not implemented for a diagonal-only band'

      outopt = 0
      outopt(1) = merge(1, 0, opt%trace)
      outopt(2) = merge(1, 0, opt%check_jacobian)
      outopt(3) = merge(1, 0, opt%return_jacobian)

      allocate(x(n), xp(n), fp(n), gp(n), scalex(n))
      x = x0
      if (allocated(opt%scalex)) then
         if (size(opt%scalex) /= n) error stop 'solve_nleqslv: scalex has wrong length'
         scalex = opt%scalex
      else
         scalex = 1.0_dp
      end if
      if (any(.not. ieee_is_finite(scalex))) error stop 'solve_nleqslv: scalex contains non-finite values'

      call liqsiz(n, qrwsiz)
      if (qrwsiz <= 0) error stop 'solve_nleqslv: LAPACK QR workspace query failed'
      lrwork = 9*n
      ldr = n
      lrjac = merge(2*n, n, method == NLEQ_BROYDEN)
      allocate(rwork(lrwork), rcdwrk(3*n), icdwrk(n), qrwork(qrwsiz), rjac(ldr,lrjac))

      active_fn => fn
      if (present(jac)) then
         active_jac => jac
      else
         nullify(active_jac)
      end if
      solver_active = .true.

      call nwnleq(x,n,scalex,maxit,jacflg,xtol,ftol,btol,cndtol,method,global,xscalm,stepmx,delta,sigma, &
         rjac,ldr,rwork,lrwork,rcdwrk,icdwrk,qrwork,qrwsiz,nleq_fjac_bridge,nleq_fvec_bridge,outopt, &
         xp,fp,gp,njcnt,nfcnt,iter,termcd)

      solver_active = .false.
      nullify(active_fn, active_jac)

      allocate(result%x(n), result%fvec(n), result%scalex(n))
      result%x = xp; result%fvec = fp; result%scalex = scalex
      result%termcd = termcd; result%nfcnt = nfcnt; result%njcnt = njcnt; result%iter = iter
      result%message = termination_message(termcd)
      result%jac_rcond = trace_last_rcond()
      if (opt%return_jacobian) then
         allocate(result%jac(n,n))
         result%jac = rjac(:,1:n)
      end if
   end subroutine solve_nleqslv

   subroutine test_nleqslv(x0, fn, result, methods, globals, nrep, options, jac)
      real(dp), intent(in) :: x0(:)
      procedure(nleq_function) :: fn
      type(nleq_test_result), intent(out) :: result
      integer, intent(in), optional :: methods(:), globals(:), nrep
      type(nleq_options), intent(in), optional :: options
      procedure(nleq_jacobian), optional :: jac
      integer, allocatable :: meth(:), glob(:)
      type(nleq_options) :: opt
      type(nleq_result) :: sol
      integer :: i, j, k, r, nr
      real(dp) :: t0, t1

      if (present(methods)) then
         allocate(meth(size(methods))); meth=methods
      else
         meth=[NLEQ_NEWTON,NLEQ_BROYDEN]
      end if
      if (present(globals)) then
         allocate(glob(size(globals))); glob=globals
      else
         glob=[NLEQ_CLINE,NLEQ_QLINE,NLEQ_GLINE,NLEQ_PWLDOG,NLEQ_DBLDOG,NLEQ_HOOK,NLEQ_NONE]
      end if
      nr=0; if(present(nrep)) nr=max(0,nrep)
      k=size(meth)*size(glob)
      allocate(result%method(k),result%global(k),result%termcd(k),result%nfcnt(k),result%njcnt(k), &
               result%iter(k),result%fnorm(k),result%cpu_seconds(k))
      k=0
      do i=1,size(meth)
         do j=1,size(glob)
            k=k+1; opt=nleq_options(); if(present(options)) opt=options
            opt%method=meth(i); opt%global=glob(j)
            call cpu_time(t0)
            if(nr>=1) then
               do r=1,nr
                  if(present(jac)) then
                     call solve_nleqslv(x0,fn,sol,opt,jac)
                  else
                     call solve_nleqslv(x0,fn,sol,opt)
                  end if
               end do
            else
               if(present(jac)) then
                  call solve_nleqslv(x0,fn,sol,opt,jac)
               else
                  call solve_nleqslv(x0,fn,sol,opt)
               end if
            end if
            call cpu_time(t1)
            result%method(k)=meth(i); result%global(k)=glob(j); result%termcd(k)=sol%termcd
            result%nfcnt(k)=sol%nfcnt; result%njcnt(k)=sol%njcnt; result%iter(k)=sol%iter
            result%fnorm(k)=0.5_dp*dot_product(sol%fvec,sol%fvec)
            result%cpu_seconds(k)=merge(t1-t0,0.0_dp,nr>=1)
         end do
      end do
   end subroutine test_nleqslv

   subroutine search_zeros(starts, fn, result, digits, options, jac)
      real(dp), intent(in) :: starts(:,:)
      procedure(nleq_function) :: fn
      type(search_zeros_result), intent(out) :: result
      integer, intent(in), optional :: digits
      type(nleq_options), intent(in), optional :: options
      procedure(nleq_jacobian), optional :: jac
      type(nleq_result) :: sol
      type(nleq_options) :: opt
      integer :: nr, nc, i, j, ndig, ncvg, nxtol, nfail, nuniq
      real(dp), allocatable :: allx(:,:), fnorm(:), rounded(:,:)
      integer, allocatable :: cvg(:), xt(:), fail(:), keep(:)
      real(dp) :: scale
      logical :: duplicate

      nr = size(starts,1); nc = size(starts,2)
      if (nr < 1 .or. nc < 1) error stop 'search_zeros: starts must be nonempty'
      ndig = 4; if (present(digits)) ndig = digits
      scale = 10.0_dp**ndig
      opt = nleq_options(); if (present(options)) opt = options
      allocate(allx(nr,nc), fnorm(nr), rounded(nr,nc), cvg(nr), xt(nr), fail(nr), keep(nr))
      allx = 0.0_dp; fnorm = huge(1.0_dp); ncvg=0; nxtol=0; nfail=0

      do i=1,nr
         if (present(jac)) then
            call solve_nleqslv(starts(i,:), fn, sol, opt, jac)
         else
            call solve_nleqslv(starts(i,:), fn, sol, opt)
         end if
         allx(i,:) = sol%x
         fnorm(i) = 0.5_dp*dot_product(sol%fvec,sol%fvec)
         if (sol%termcd == 1) then
            ncvg=ncvg+1; cvg(ncvg)=i
         else if (sol%termcd == 2) then
            nxtol=nxtol+1; xt(nxtol)=i
         else
            nfail=nfail+1; fail(nfail)=i
         end if
      end do

      allocate(result%converged_rows(ncvg), result%xtol_rows(nxtol), result%failed_rows(nfail))
      if(ncvg>0) result%converged_rows=cvg(:ncvg)
      if(nxtol>0) result%xtol_rows=xt(:nxtol)
      if(nfail>0) result%failed_rows=fail(:nfail)
      if (ncvg == 0) then
         allocate(result%x(0,nc), result%xfnorm(0)); return
      end if

      do i=1,ncvg
         rounded(i,:) = anint(allx(cvg(i),:)*scale)/scale
      end do
      nuniq=0
      do i=1,ncvg
         duplicate=.false.
         do j=1,nuniq
            if (all(rounded(i,:) == rounded(keep(j),:))) then
               duplicate=.true.; exit
            end if
         end do
         if(.not.duplicate) then
            nuniq=nuniq+1; keep(nuniq)=i
         end if
      end do
      allocate(result%x(nuniq,nc), result%xfnorm(nuniq))
      do j=1,nuniq
         i=cvg(keep(j)); result%x(j,:)=allx(i,:); result%xfnorm(j)=fnorm(i)
      end do
      call sort_solutions(result%x,result%xfnorm)
   end subroutine search_zeros

   subroutine sort_solutions(x, fnorm)
      real(dp), intent(inout) :: x(:,:), fnorm(:)
      integer :: i,j,k,nc
      real(dp), allocatable :: tmp(:)
      real(dp) :: tf
      logical :: less
      nc=size(x,2); allocate(tmp(nc))
      do i=2,size(x,1)
         tmp=x(i,:); tf=fnorm(i); j=i-1
         do while(j>=1)
            less=.false.
            do k=1,nc
               if(tmp(k)<x(j,k)) then; less=.true.; exit
               else if(tmp(k)>x(j,k)) then; exit
               end if
            end do
            if(.not.less) exit
            x(j+1,:)=x(j,:); fnorm(j+1)=fnorm(j); j=j-1
         end do
         x(j+1,:)=tmp; fnorm(j+1)=tf
      end do
   end subroutine sort_solutions

   subroutine nleq_fvec_bridge(x, f, n, flag)
      integer, intent(in) :: n, flag
      real(dp), intent(in) :: x(*)
      real(dp), intent(out) :: f(*)
      real(dp) :: xt(n), ft(n)
      integer :: i
      if (.not.associated(active_fn)) error stop 'nleqslv internal: no active function callback'
      xt=x(1:n); call active_fn(xt,ft); f(1:n)=ft
      if (any(.not.ieee_is_finite(ft))) then
         if (flag /= 0) error stop 'solve_nleqslv: non-finite function value during numerical Jacobian'
         do concurrent(i=1:n)
            if(.not.ieee_is_finite(f(i))) f(i)=sqrt(huge(1.0_dp)/real(n,dp))
         end do
      end if
   end subroutine nleq_fvec_bridge

   subroutine nleq_fjac_bridge(rjac, ldr, x, n)
      integer, intent(in) :: ldr, n
      real(dp), intent(out) :: rjac(ldr,*)
      real(dp), intent(in) :: x(*)
      real(dp) :: xt(n), jt(n,n)
      if (.not.associated(active_jac)) error stop 'nleqslv internal: no active Jacobian callback'
      xt=x(1:n); call active_jac(xt,jt)
      if(any(.not.ieee_is_finite(jt))) error stop 'solve_nleqslv: non-finite user Jacobian'
      rjac(1:n,1:n)=jt
   end subroutine nleq_fjac_bridge

   function termination_message(termcd) result(msg)
      integer, intent(in) :: termcd
      character(len=160) :: msg
      select case(termcd)
      case(1); msg='Function criterion near zero'
      case(2); msg="x-values within tolerance 'xtol'"
      case(3); msg='No better point found (algorithm has stalled)'
      case(4); msg='Iteration limit exceeded'
      case(5); msg='Jacobian is too ill-conditioned'
      case(6); msg='Jacobian is singular'
      case(7); msg='Jacobian is completely unusable'
      case(-10); msg='User supplied Jacobian most likely incorrect'
      case default; write(msg,'("Termination code ",i0)') termcd
      end select
   end function termination_message


   real(dp) function trace_last_rcond()
      use nleqslv_trace_state, only : last_rcond
      trace_last_rcond = last_rcond
   end function trace_last_rcond

end module nleqslv_fortran
