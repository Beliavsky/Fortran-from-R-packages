! SPDX-License-Identifier: GPL-2.0-only
module nlsic_solver
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use nlsic_kinds, only : dp
   use nlsic_types
   use nlsic_linalg, only : pseudoinverse, null_space, vecnorm, student_t_quantile
   use nlsic_nnls, only : ldp_solve
   use nlsic_linear, only : lsi, lsi_ln, lsie_ln
   implicit none
   private
   public :: nlsic_solve, numerical_jacobian
contains

   subroutine numerical_jacobian(par,m,residual,jac,ierr)
      real(dp), intent(in) :: par(:)
      integer, intent(in) :: m
      procedure(residual_function) :: residual
      real(dp), intent(out) :: jac(:,:)
      integer, intent(out) :: ierr
      real(dp), allocatable :: xp(:),xm(:),rp(:),rm(:)
      real(dp) :: h
      integer :: n,j,ip,im
      n=size(par); ierr=0; jac=0.0_dp
      if(size(jac,1)/=m .or. size(jac,2)/=n) then; ierr=-1; return; end if
      allocate(xp(n),xm(n),rp(m),rm(m))
      do j=1,n
         h=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(par(j)))
         xp=par; xm=par; xp(j)=xp(j)+h; xm(j)=xm(j)-h
         call residual(xp,rp,ip); call residual(xm,rm,im)
         if(ip/=0 .or. im/=0) then; ierr=1; return; end if
         jac(:,j)=(rp-rm)/(2.0_dp*h)
      end do
   end subroutine numerical_jacobian

   subroutine eval_residual(par,m,residual,res,ierr)
      real(dp), intent(in) :: par(:)
      integer, intent(in) :: m
      procedure(residual_function) :: residual
      real(dp), intent(out) :: res(:)
      integer, intent(out) :: ierr
      if(size(res)/=m) then; ierr=-1; return; end if
      call residual(par,res,ierr)
   end subroutine eval_residual

   subroutine eval_jacobian(par,m,residual,res,jac,ierr,jacobian)
      real(dp), intent(in) :: par(:)
      integer, intent(in) :: m
      procedure(residual_function) :: residual
      real(dp), intent(out) :: res(:),jac(:,:)
      integer, intent(out) :: ierr
      procedure(jacobian_function), optional :: jacobian
      if(present(jacobian)) then
         call jacobian(par,res,jac,ierr)
      else
         call residual(par,res,ierr)
         if(ierr/=0) return
         call numerical_jacobian(par,m,residual,jac,ierr)
      end if
   end subroutine eval_jacobian

   subroutine nlsic_solve(par0,m,residual,result,jacobian,u,co,e,eco,control,mnorm)
      real(dp), intent(in) :: par0(:)
      integer, intent(in) :: m
      procedure(residual_function) :: residual
      type(nlsic_result), intent(out) :: result
      procedure(jacobian_function), optional :: jacobian
      real(dp), intent(in), optional :: u(:,:),co(:),e(:,:),eco(:),mnorm(:,:)
      type(nlsic_control), intent(in), optional :: control
      type(nlsic_control) :: con
      type(lsi_result) :: lr
      real(dp), allocatable :: par(:),r(:),rold(:),jmat(:,:),jvalid(:,:),rv(:),jstep(:)
      real(dp), allocatable :: u0(:,:),c0(:),e0(:,:),ec0(:),ine(:),step(:),laststep(:)
      real(dp), allocatable :: ep(:,:),bn(:,:),dz(:),pinvj(:,:),covz(:,:),jz(:,:),zeros_e(:)
      logical, allocatable :: valid(:)
      real(dp) :: tol,norb,nornew,normp,k,klast,kms,actual,pred,resap,n2ap,hfac,tcrit
      real(dp) :: fdeg,scale_co
      integer :: n,nu,ne,it,btit,reuse_count,st,ranke,rankj,nvalid,slot,ierr
      logical :: converged,feasible,descending,force_jac,accepted
      character(:), allocatable :: final_message

      n=size(par0); con=nlsic_control(); if(present(control)) con=control
      call allocate_result_shell(result,n,m,con)
      if(n<=0 .or. m<=0 .or. con%maxit<1 .or. con%btmaxit<1 .or. con%rcond<=1.0_dp) then
         call set_error(result,NLSIC_INVALID_INPUT,'nlsic: invalid dimensions or control values'); return
      end if
      if(con%btfrac<=0.0_dp .or. con%btfrac>=1.0_dp .or. con%btdesc<0.0_dp .or. con%btdesc>1.0_dp) then
         call set_error(result,NLSIC_INVALID_INPUT,'nlsic: invalid backtracking control'); return
      end if
      if(con%maxstep>0.0_dp .and. con%maxstep<=con%errx) then
         call set_error(result,NLSIC_INVALID_INPUT,'nlsic: maxstep must exceed errx'); return
      end if
      tol=1.0_dp/con%rcond
      call copy_constraints(n,u,co,u0,c0,st)
      if(st/=0) then; call set_error(result,NLSIC_INVALID_INPUT,'nlsic: bad inequality dimensions'); return; end if
      call copy_equalities(n,e,eco,e0,ec0,st)
      if(st/=0) then; call set_error(result,NLSIC_INVALID_INPUT,'nlsic: bad equality dimensions'); return; end if
      nu=size(u0,1); ne=size(e0,1)
      allocate(par(n),r(m),rold(m),jmat(m,n),valid(m),step(n),laststep(n))
      par=par0; laststep=0.0_dp; step=0.0_dp

      ! Project initial point to the equality manifold, preserving its null-space component.
      if(ne>0) then
         if(ne>n) then; call set_error(result,NLSIC_INVALID_INPUT,'nlsic: equality system is overdetermined'); return; end if
         allocate(ep(n,ne)); call pseudoinverse(e0,ep,ranke,con%rcond,st)
         if(ranke<ne) then; call set_error(result,NLSIC_INVALID_INPUT,'nlsic: equality matrix is rank deficient'); return; end if
         call null_space(e0,bn,ranke,con%rcond,st)
         par=matmul(ep,ec0)
         if(size(bn,2)>0) par=par+matmul(bn,matmul(transpose(bn),par0-par))
         if(maxval(abs(matmul(e0,par)-ec0))>100.0_dp*tol*max(1.0_dp,maxval(abs(ec0)))) then
            call set_error(result,NLSIC_INFEASIBLE,'nlsic: inconsistent equality constraints'); return
         end if
      else
         allocate(bn(n,n)); bn=0.0_dp
         do st=1,n; bn(st,st)=1.0_dp; end do
      end if

      if(nu>0) then
         allocate(ine(nu)); ine=c0-matmul(u0,par)
         if(any(ine>tol)) then
            if(size(bn,2)==0) then
               call set_error(result,NLSIC_INFEASIBLE,'nlsic: infeasible starting point'); return
            end if
            allocate(dz(size(bn,2)))
            call ldp_solve(matmul(u0,bn),ine,dz,feasible,con%rcond)
            if(.not.feasible) then
               call set_error(result,NLSIC_INFEASIBLE,'nlsic: infeasible constraints at starting point'); return
            end if
            par=par+matmul(bn,dz); deallocate(dz)
         end if
      end if

      call eval_jacobian(par,m,residual,r,jmat,ierr,jacobian)
      if(ierr/=0) then; call set_error(result,NLSIC_NUMERICAL,'nlsic: initial residual/Jacobian evaluation failed'); return; end if
      valid=ieee_is_finite(r); nvalid=count(valid)
      if(nvalid==0) then; call set_error(result,NLSIC_NUMERICAL,'nlsic: no finite residual value'); return; end if
      allocate(rv(nvalid),jvalid(nvalid,n)); rv=pack(r,valid); jvalid=reshape(pack(jmat,spread(valid,2,n)),[nvalid,n])
      norb=vecnorm(rv); rold=r
      if(con%history) then
         result%par_history(:,1)=par; result%residual_history(:,1)=r
      end if
      converged=norb<1.0e-20_dp; it=0; btit=0; reuse_count=0; force_jac=.false.; normp=0.0_dp
      final_message=''

      do while(.not.converged .and. it<con%maxit)
         if(it>0) then
            if(force_jac .or. .not.con%reuse_jac .or. btit>1 .or. reuse_count>=con%max_reuse) then
               call eval_jacobian(par,m,residual,r,jmat,ierr,jacobian)
               if(ierr/=0) then; call set_error(result,NLSIC_NUMERICAL,'nlsic: residual/Jacobian evaluation failed'); return; end if
               valid=ieee_is_finite(r); nvalid=count(valid)
               if(nvalid==0) then; call set_error(result,NLSIC_NUMERICAL,'nlsic: no finite residual value'); return; end if
               deallocate(rv,jvalid); allocate(rv(nvalid),jvalid(nvalid,n))
               rv=pack(r,valid); jvalid=reshape(pack(jmat,spread(valid,2,n)),[nvalid,n])
               reuse_count=0; force_jac=.false.; norb=vecnorm(rv)
            else
               reuse_count=reuse_count+1
            end if
         end if
         rold=r
         if(nu>0) then
            if(allocated(ine)) deallocate(ine); allocate(ine(nu)); ine=c0-matmul(u0,par)
         end if
         if(ne>0) then
            allocate(zeros_e(ne)); zeros_e=0.0_dp
            if(con%least_norm_step) then
               if(present(mnorm)) then
                  call lsie_ln(jvalid,-rv,lr,u0,ine,e0,zeros_e,con%rcond,mnorm,-par)
               else
                  call lsie_ln(jvalid,-rv,lr,u0,ine,e0,zeros_e,con%rcond,x0=-par)
               end if
            else
               call lsie_ln(jvalid,-rv,lr,u0,ine,e0,zeros_e,con%rcond)
            end if
            deallocate(zeros_e)
         else if(con%least_norm_step) then
            if(present(mnorm)) then
               call lsi_ln(jvalid,-rv,lr,u0,ine,con%rcond,mnorm,-par)
            else
               call lsi_ln(jvalid,-rv,lr,u0,ine,con%rcond,x0=-par)
            end if
         else
            call lsi(jvalid,-rv,lr,u0,ine,con%rcond)
         end if
         if(.not.lr%succeeded()) then
            call set_error(result,NLSIC_NUMERICAL,'nlsic: linearized constrained least-squares step failed'); return
         end if
         step=lr%x; normp=vecnorm(step)
         if(nu>0) then
            scale_co=max(1.0_dp,maxval(abs(c0)))
            if(minval(matmul(u0,par+step)-c0)<-1000.0_dp*tol*scale_co) then
               call set_error(result,NLSIC_NUMERICAL,'nlsic: linearized step violates inequalities'); return
            end if
         end if
         if(ne>0) then
            if(maxval(abs(matmul(e0,step)))>1000.0_dp*tol) then
               call set_error(result,NLSIC_NUMERICAL,'nlsic: linearized step violates equalities'); return
            end if
         end if
         if(con%history) then
            slot=min(it+1,size(result%direction_history,2)); result%direction_history(:,slot)=step
            result%reuse_history(slot)=reuse_count
         end if
         converged=normp<=con%errx
         if(normp<=tiny(1.0_dp)) then
            laststep=0.0_dp; it=it+1; btit=0
            call eval_residual(par,m,residual,r,ierr)
            if(ierr/=0) then; call set_error(result,NLSIC_NUMERICAL,'nlsic: residual evaluation failed'); return; end if
            exit
         end if
         allocate(jstep(nvalid)); jstep=matmul(jvalid,step); resap=dot_product(rv,jstep); n2ap=dot_product(jstep,jstep)
         if(resap>0.0_dp) then
            deallocate(jstep)
            call set_error(result,NLSIC_NOT_DESCENT, &
               'nlsic: linearized solver returned a non-descending direction')
            return
         end if
         k=con%btstart
         if(con%maxstep>0.0_dp) then; kms=con%maxstep/normp; k=min(k,kms); end if
         btit=0; descending=.false.; accepted=.false.; klast=k
         do while(.not.descending .and. btit<con%btmaxit .and. k>=con%btkmin)
            laststep=k*step
            call eval_residual(par+laststep,m,residual,r,ierr)
            if(ierr/=0) then
               k=max(k*con%btfrac,con%btkmin); btit=btit+1; cycle
            end if
            if(any(.not.ieee_is_finite(pack(r,valid)))) then
               k=max(k*con%btfrac,con%btkmin); btit=btit+1; cycle
            end if
            nornew=vecnorm(pack(r,valid))
            actual=norb*norb-nornew*nornew
            pred=-k*(2.0_dp*resap+k*n2ap)
            descending=actual>=con%btdesc*pred
            klast=k
            if(con%adaptbt .and. abs(resap*k)>tiny(1.0_dp)) then
               hfac=max(min(2.0_dp,1.0_dp/(2.0_dp-(nornew-norb)*(nornew+norb)/(resap*k))),con%btfrac)
               if(hfac>1.0_dp .and. nu>0) then
                  if(any(matmul(u0,par+(k*hfac)*step)-c0 < -tol)) hfac=con%btfrac
               end if
               k=max(k*hfac,con%btkmin)
            else
               k=max(k*con%btfrac,con%btkmin)
            end if
            btit=btit+1
            if(descending) accepted=.true.
         end do
         deallocate(jstep)
         if(.not.accepted) then
            if(reuse_count>0) then
               force_jac=.true.; call eval_residual(par,m,residual,r,ierr)
               if(ierr/=0) then; call set_error(result,NLSIC_NUMERICAL,'nlsic: residual restart failed'); return; end if
               cycle
            end if
            call set_error(result,NLSIC_BACKTRACK_LIMIT,'nlsic: backtracking failed to find an acceptable step'); return
         end if
         par=par+laststep; it=it+1
         if(con%monotone .and. nornew>norb) then
            converged=.true.; final_message='nlsic: non-monotone cost increase encountered'
         end if
         norb=nornew
         if(con%history) then
            result%step_history(:,it)=laststep
            result%par_history(:,it+1)=par; result%residual_history(:,it+1)=r
         end if
      end do

      result%par=par; result%lastp=step; result%laststep=laststep; result%normp=normp
      result%residuals=r; result%previous_residuals=rold; result%jacobian=jmat
      result%iterations=it; result%backtrack_iterations=btit; result%ci_p=con%ci_p
      result%converged=converged
      if(it>=con%maxit .and. .not.converged) then
         result%status=NLSIC_MAX_ITER; result%message='nlsic: maximal nonlinear iteration number reached'
      else
         result%status=NLSIC_SUCCESS; if(len(final_message)>0) result%message=final_message
      end if

      if(con%report_ci .and. result%status==NLSIC_SUCCESS) then
         call eval_jacobian(par,m,residual,r,jmat,ierr,jacobian)
         if(ierr==0) then
            valid=ieee_is_finite(r); nvalid=count(valid)
            if(ne>0) then
               jz=matmul(reshape(pack(jmat,spread(valid,2,n)),[nvalid,n]),bn)
               allocate(pinvj(size(jz,2),size(jz,1))); call pseudoinverse(jz,pinvj,rankj,con%rcond,st)
               covz=matmul(pinvj,transpose(pinvj)); result%covariance=matmul(bn,matmul(covz,transpose(bn)))
            else
               jvalid=reshape(pack(jmat,spread(valid,2,n)),[nvalid,n])
               allocate(pinvj(n,nvalid)); call pseudoinverse(jvalid,pinvj,rankj,con%rcond,st)
               result%covariance=matmul(pinvj,transpose(pinvj))
            end if
            fdeg=max(1.0_dp,real(nvalid-n+ne,dp)); result%ci_fdeg=fdeg
            result%sd_res=vecnorm(pack(r,valid))/sqrt(fdeg)
            tcrit=student_t_quantile(0.5_dp*(1.0_dp+con%ci_p),fdeg)
            do st=1,n
               result%hci(st)=tcrit*result%sd_res*sqrt(max(0.0_dp,result%covariance(st,st)))
            end do
            result%jacobian=jmat; result%residuals=r
         end if
      end if
      if(con%history) call trim_history(result,it,m,n)
   end subroutine nlsic_solve

   subroutine allocate_result_shell(result,n,m,con)
      type(nlsic_result), intent(out) :: result
      integer, intent(in) :: n,m
      type(nlsic_control), intent(in) :: con
      allocate(result%par(n),result%lastp(n),result%laststep(n),result%residuals(m), &
         result%previous_residuals(m),result%jacobian(m,n),result%covariance(n,n),result%hci(n))
      result%par=0.0_dp; result%lastp=0.0_dp; result%laststep=0.0_dp; result%residuals=0.0_dp
      result%previous_residuals=0.0_dp; result%jacobian=0.0_dp; result%covariance=0.0_dp; result%hci=0.0_dp
      if(con%history) then
         allocate(result%par_history(n,con%maxit+1),result%direction_history(n,con%maxit), &
            result%step_history(n,con%maxit),result%residual_history(m,con%maxit+1), &
            result%reuse_history(con%maxit))
         result%par_history=0.0_dp; result%direction_history=0.0_dp; result%step_history=0.0_dp
         result%residual_history=0.0_dp; result%reuse_history=0
      end if
   end subroutine allocate_result_shell

   subroutine set_error(result,status,message)
      type(nlsic_result), intent(inout) :: result
      integer, intent(in) :: status
      character(*), intent(in) :: message
      result%status=status; result%message=message
   end subroutine set_error

   subroutine copy_constraints(n,u,co,u0,c0,status)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: u(:,:),co(:)
      real(dp), allocatable, intent(out) :: u0(:,:),c0(:)
      integer, intent(out) :: status
      status=0
      if(present(u)) then
         if(size(u,2)/=n .or. .not.present(co)) then; status=1; allocate(u0(0,n),c0(0)); return; end if
         if(size(co)/=size(u,1)) then; status=1; allocate(u0(0,n),c0(0)); return; end if
         allocate(u0(size(u,1),n),c0(size(co))); u0=u; c0=co
      else
         allocate(u0(0,n),c0(0))
      end if
   end subroutine copy_constraints

   subroutine copy_equalities(n,e,eco,e0,ec0,status)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: e(:,:),eco(:)
      real(dp), allocatable, intent(out) :: e0(:,:),ec0(:)
      integer, intent(out) :: status
      status=0
      if(present(e)) then
         if(size(e,2)/=n .or. .not.present(eco)) then; status=1; allocate(e0(0,n),ec0(0)); return; end if
         if(size(eco)/=size(e,1)) then; status=1; allocate(e0(0,n),ec0(0)); return; end if
         allocate(e0(size(e,1),n),ec0(size(eco))); e0=e; ec0=eco
      else
         allocate(e0(0,n),ec0(0))
      end if
   end subroutine copy_equalities

   subroutine trim_history(result,it,m,n)
      type(nlsic_result), intent(inout) :: result
      integer, intent(in) :: it,m,n
      real(dp), allocatable :: ph(:,:),dh(:,:),sh(:,:),rh(:,:)
      integer, allocatable :: uh(:)
      allocate(ph(n,it+1),rh(m,it+1)); ph=result%par_history(:,1:it+1); rh=result%residual_history(:,1:it+1)
      call move_alloc(ph,result%par_history); call move_alloc(rh,result%residual_history)
      allocate(dh(n,it),sh(n,it),uh(it))
      if(it>0) then
         dh=result%direction_history(:,1:it); sh=result%step_history(:,1:it); uh=result%reuse_history(1:it)
      end if
      call move_alloc(dh,result%direction_history)
      call move_alloc(sh,result%step_history)
      call move_alloc(uh,result%reuse_history)
   end subroutine trim_history
end module nlsic_solver
