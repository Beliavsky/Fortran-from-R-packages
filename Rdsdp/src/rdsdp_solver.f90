! Dual-scaling/log-barrier SDP solver translated from DSDP5 concepts.
!
! v0.2.0 added abstract dense/sparse/low-rank SDP data and sparse Schur assembly.
! v0.3.0 adds matrix-free PCG plus RCM-ordered sparse LDL^T Schur solves.
! License/provenance: licenses/DSDP-LICENSE.
module rdsdp_solver
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rdsdp_kinds, only : dp
   use rdsdp_types, only : dsdp_problem, dsdp_control, dsdp_solution, &
      dsdp_sdp_block, dsdp_lp_block, dsdp_pdfeasible, dsdp_unknown, dsdp_infeasible, &
      dsdp_converged, dsdp_max_iterations, dsdp_linesearch_failure, dsdp_singular_schur, &
      dsdp_numeric_failure, dsdp_data_dense, dsdp_data_sparse, dsdp_data_lowrank
   use rdsdp_linalg, only : spd_inverse_logdet, spd_logdet, solve_spd, norm2_dp, min_eigenvalue_sym
   use rdsdp_problem_mod, only : validate_problem, compress_problem, densify_problem
   use rdsdp_data, only : data_add_scaled, data_dot, data_fnorm2, data_nnz, data_schur_pair, &
      data_storage, get_data_dense
   use rdsdp_sparse_ldlt, only : sparse_ldlt_cache, sparse_ldlt_solve
   implicit none
   private
   public :: dsdp_solve

   type :: newton_block_context
      integer :: category = dsdp_sdp_block
      integer :: n = 0
      real(dp), allocatable :: sinv(:,:)
      real(dp), allocatable :: w2(:)
   end type newton_block_context

   type :: newton_context
      integer :: m = 0
      integer :: nv = 0
      real(dp) :: mu = 0.0_dp
      real(dp) :: r = 0.0_dp
      real(dp), allocatable :: diag(:)
      type(newton_block_context), allocatable :: block(:)
   end type newton_context

contains

   subroutine dsdp_solve(prob,sol,control,y0)
      type(dsdp_problem), intent(in) :: prob
      type(dsdp_solution), intent(out) :: sol
      type(dsdp_control), intent(in), optional :: control
      real(dp), intent(in), optional :: y0(:)

      type(dsdp_control) :: ctrl
      type(dsdp_problem) :: work
      type(newton_context) :: nctx
      type(sparse_ldlt_cache) :: scache
      real(dp), allocatable :: y(:),z(:),g(:),h(:,:),dz(:),zt(:)
      real(dp) :: r,mu,phi,phit,dirder,alpha,pnorm
      real(dp) :: normc,normb,rscale,mu_floor,t0,t1
      integer :: outer,inner,ls,nv,nu,status
      logical :: ok,accepted,use_mf
      character(len=:), allocatable :: message

      ctrl=dsdp_control()
      if (present(control)) ctrl=control
      call validate_problem(prob,ok,message)
      if (.not.ok) error stop 'dsdp_solve: '//message
      if (ctrl%penalty<=0.0_dp) error stop 'dsdp_solve: penalty must be positive'
      if (ctrl%mu_factor<=0.0_dp .or. ctrl%mu_factor>=1.0_dp) &
         error stop 'dsdp_solve: mu_factor must be in (0,1)'
      if (ctrl%cg_maxiter<1) error stop 'dsdp_solve: cg_maxiter must be positive'

      work=prob
      if (ctrl%use_sparse_data) then
         call compress_problem(work,ctrl%sparse_density_threshold,.true.)
      else
         call densify_problem(work)
      end if

      sol=dsdp_solution()
      sol%sdp_data_nnz=count_problem_nnz(work)
      nv=work%m+1
      nu=work%barrier_dimension()
      allocate(y(work%m),z(nv),g(nv),h(nv,nv),dz(nv),zt(nv))
      y=0.0_dp
      if (present(y0)) then
         if (size(y0)/=work%m) error stop 'dsdp_solve: y0 has wrong size'
         y=y0
      end if
      normc=data_c_norm(work)
      normb=norm2_dp(work%b)
      r=initial_shift(work,y)
      z(1:work%m)=y; z(nv)=r
      mu=max(1.0_dp,ctrl%penalty*r/max(1.0_dp,real(nu,dp)))
      mu=max(mu,0.1_dp*(1.0_dp+normb+normc))
      mu_floor=max(1.0e-14_dp,ctrl%gaptol*1.0e-3_dp)
      rscale=1.0_dp+normc
      status=dsdp_max_iterations

      do outer=1,ctrl%maxiter
         do inner=1,ctrl%max_correctors
            use_mf=ctrl%use_cg .and. ctrl%cg_matrix_free .and. nv>=max(1,ctrl%cg_threshold)
            call cpu_time(t0)
            if (use_mf) then
               call evaluate_newton_matrixfree(work,z,mu,ctrl%penalty,phi,g,nctx,sol,ok)
            else
               call evaluate_newton(work,z,mu,ctrl%penalty,phi,g,h,sol,ok)
            end if
            call cpu_time(t1); sol%schur_assembly_time=sol%schur_assembly_time+(t1-t0)
            if (.not.ok .or. .not.ieee_is_finite(phi)) then
               status=dsdp_numeric_failure; exit
            end if
            call cpu_time(t0)
            if (use_mf) then
               call pcg_matrixfree(work,nctx,g,dz,ctrl,sol,ok)
               if (.not.ok .and. ctrl%cg_fallback_direct) then
                  call evaluate_newton(work,z,mu,ctrl%penalty,phi,g,h,sol,ok)
                  if (ok) call solve_schur_system(h,g,dz,ctrl,sol,scache,ok)
               end if
            else
               call solve_schur_system(h,g,dz,ctrl,sol,scache,ok)
            end if
            call cpu_time(t1); sol%schur_solve_time=sol%schur_solve_time+(t1-t0)
            if (.not.ok) then
               status=dsdp_singular_schur; exit
            end if
            dirder=dot_product(g,dz)
            if (dirder<0.0_dp .or. .not.ieee_is_finite(dirder)) then
               status=dsdp_numeric_failure; exit
            end if
            pnorm=sqrt(max(0.0_dp,dirder/max(mu,tiny(1.0_dp))))
            sol%pnorm=pnorm
            if (pnorm<=ctrl%newton_tol) exit

            alpha=1.0_dp; accepted=.false.
            do ls=0,ctrl%max_linesearch-1
               zt=z+alpha*dz
               call state_value(work,zt,mu,ctrl%penalty,phit,ok)
               if (ok) then
                  if (phit>=phi+ctrl%armijo*alpha*dirder) then
                     accepted=.true.; exit
                  end if
               end if
               alpha=0.5_dp*alpha
               sol%line_search_backtracks=sol%line_search_backtracks+1
            end do
            if (.not.accepted) then
               status=dsdp_linesearch_failure; exit
            end if
            z=zt; sol%dstep=alpha; sol%newton_steps=sol%newton_steps+1
         end do
         if (status/=dsdp_max_iterations) exit

         y=z(1:work%m); r=z(nv)
         call recover_primal_and_metrics(work,y,r,mu,sol,ok)
         if (.not.ok) then
            status=dsdp_numeric_failure; exit
         end if
         sol%iterations=outer; sol%mu=mu; sol%r=r; sol%y=y

         if (ctrl%print>0) then
            write(*,'("it=",i3," mu=",es10.3," r=",es10.3," pobj=",es15.7," dobj=",es15.7, &
               &" gap=",es9.2," pinf=",es9.2," pnorm=",es9.2)') outer,mu,r,sol%pobj,sol%dobj, &
               sol%relgap,sol%pinfeas,sol%pnorm
         end if

         if (sol%relgap<=ctrl%gaptol .and. sol%pinfeas<=ctrl%pinfeastol .and. &
             r<=ctrl%rtol*rscale) then
            status=dsdp_converged; exit
         end if
         mu=max(mu_floor,mu*ctrl%mu_factor)
      end do

      if (.not.allocated(sol%y)) then
         y=z(1:work%m); r=z(nv)
         call recover_primal_and_metrics(work,y,r,mu,sol,ok)
         sol%y=y; sol%r=r; sol%mu=mu
      end if
      sol%status=status
      if (status==dsdp_converged) then
         sol%stype=dsdp_pdfeasible
      else if (sol%r>1.0e-3_dp*rscale .and. sol%mu<=100.0_dp*mu_floor) then
         sol%stype=dsdp_infeasible
      else
         sol%stype=dsdp_unknown
      end if
   end subroutine dsdp_solve

   subroutine evaluate_newton(prob,z,mu,penalty,phi,g,h,sol,ok)
      type(dsdp_problem), intent(in) :: prob
      real(dp), intent(in) :: z(:),mu,penalty
      real(dp), intent(out) :: phi,g(:),h(:,:)
      type(dsdp_solution), intent(inout) :: sol
      logical, intent(out) :: ok
      real(dp), allocatable :: s(:,:),sinv(:,:),s2(:,:),sl(:),w(:)
      real(dp) :: logdet,r,cross
      integer :: k,i,j,n,m,nv
      logical :: pd
      m=prob%m; nv=m+1; r=z(nv)
      g=0.0_dp; h=0.0_dp; sol%schur_assemblies=sol%schur_assemblies+1
      g(1:m)=prob%b; g(nv)=-penalty
      phi=dot_product(prob%b,z(1:m))-penalty*r
      if (r<=0.0_dp .or. .not.ieee_is_finite(r)) then; ok=.false.; return; end if
      phi=phi+mu*log(r); g(nv)=g(nv)+mu/r; h(nv,nv)=h(nv,nv)+mu/(r*r)

      do k=1,size(prob%block)
         n=prob%block(k)%n
         select case(prob%block(k)%category)
         case(dsdp_sdp_block)
            allocate(s(n,n)); s=0.0_dp
            call data_add_scaled(prob%block(k),0,1.0_dp,s)
            do i=1,m; call data_add_scaled(prob%block(k),i,-z(i),s); end do
            do i=1,n; s(i,i)=s(i,i)+r; end do
            call spd_inverse_logdet(s,sinv,logdet,pd)
            if (.not.pd) then; ok=.false.; return; end if
            phi=phi+mu*logdet
            do i=1,m; g(i)=g(i)-mu*data_dot(prob%block(k),i,sinv); end do
            g(nv)=g(nv)+mu*trace_mat(sinv)
            s2=matmul(sinv,sinv)
            call add_sdp_hessian(prob%block(k),m,sinv,mu,h(1:m,1:m),sol)
            do j=1,m
               cross=-mu*data_dot(prob%block(k),j,s2)
               h(j,nv)=h(j,nv)+cross; h(nv,j)=h(j,nv)
            end do
            h(nv,nv)=h(nv,nv)+mu*sum(sinv*sinv)
            deallocate(s,sinv,s2)

         case(dsdp_lp_block)
            allocate(sl(n),w(n))
            sl=prob%block(k)%cdiag-matmul(prob%block(k)%adiag,z(1:m))+r
            if (any(sl<=0.0_dp) .or. any(.not.ieee_is_finite(sl))) then; ok=.false.; return; end if
            phi=phi+mu*sum(log(sl)); w=1.0_dp/sl
            do i=1,m; g(i)=g(i)-mu*dot_product(prob%block(k)%adiag(:,i),w); end do
            g(nv)=g(nv)+mu*sum(w); w=w*w
            do j=1,m
               do i=1,j
                  h(i,j)=h(i,j)+mu*dot_product(prob%block(k)%adiag(:,i)*w,prob%block(k)%adiag(:,j))
                  if (i/=j) h(j,i)=h(i,j)
               end do
               cross=-mu*dot_product(prob%block(k)%adiag(:,j),w)
               h(j,nv)=h(j,nv)+cross; h(nv,j)=h(j,nv)
            end do
            h(nv,nv)=h(nv,nv)+mu*sum(w); deallocate(sl,w)
         end select
      end do
      ok=all(ieee_is_finite(g)) .and. all(ieee_is_finite(h)) .and. ieee_is_finite(phi)
   end subroutine evaluate_newton

   subroutine add_sdp_hessian(block,m,sinv,mu,h,sol)
      use rdsdp_types, only : dsdp_block
      type(dsdp_block), intent(in) :: block
      integer, intent(in) :: m
      real(dp), intent(in) :: sinv(:,:),mu
      real(dp), intent(inout) :: h(:,:)
      type(dsdp_solution), intent(inout) :: sol
      real(dp), allocatable :: bdense(:,:,:),a(:,:)
      real(dp) :: v
      integer :: i,j,n,si,sj
      logical :: anydense
      n=block%n; anydense=.false.
      do i=1,m
         if (data_storage(block,i)==dsdp_data_dense) anydense=.true.
      end do
      if (anydense) then
         allocate(bdense(n,n,m)); bdense=0.0_dp
         do i=1,m
            if (data_storage(block,i)==dsdp_data_dense) then
               call get_data_dense(block,i,a)
               bdense(:,:,i)=matmul(matmul(sinv,a),sinv)
            end if
         end do
         do j=1,m
            sj=data_storage(block,j)
            do i=1,j
               si=data_storage(block,i)
               if (si==dsdp_data_dense) then
                  v=data_dot(block,j,bdense(:,:,i))
               else if (sj==dsdp_data_dense) then
                  v=data_dot(block,i,bdense(:,:,j))
               else
                  v=data_schur_pair(block,i,j,sinv)
               end if
               if (si==dsdp_data_dense .or. sj==dsdp_data_dense) then
                  sol%dense_pair_evals=sol%dense_pair_evals+1
               else if (si==dsdp_data_lowrank .or. sj==dsdp_data_lowrank) then
                  sol%lowrank_pair_evals=sol%lowrank_pair_evals+1
               else
                  sol%sparse_pair_evals=sol%sparse_pair_evals+1
               end if
               h(i,j)=h(i,j)+mu*v
               if (i/=j) h(j,i)=h(i,j)
            end do
         end do
      else
         do j=1,m
            sj=data_storage(block,j)
            do i=1,j
               si=data_storage(block,i)
               v=data_schur_pair(block,i,j,sinv)
               if (si==dsdp_data_lowrank .or. sj==dsdp_data_lowrank) then
                  sol%lowrank_pair_evals=sol%lowrank_pair_evals+1
               else
                  sol%sparse_pair_evals=sol%sparse_pair_evals+1
               end if
               h(i,j)=h(i,j)+mu*v
               if (i/=j) h(j,i)=h(i,j)
            end do
         end do
      end if
   end subroutine add_sdp_hessian

   subroutine evaluate_newton_matrixfree(prob,z,mu,penalty,phi,g,ctx,sol,ok)
      type(dsdp_problem), intent(in) :: prob
      real(dp), intent(in) :: z(:),mu,penalty
      real(dp), intent(out) :: phi,g(:)
      type(newton_context), intent(out) :: ctx
      type(dsdp_solution), intent(inout) :: sol
      logical, intent(out) :: ok
      real(dp), allocatable :: s(:,:),sinv(:,:),sl(:),w(:)
      real(dp) :: logdet,r,v
      integer :: k,i,n,m,nv,storage
      logical :: pd

      m=prob%m; nv=m+1; r=z(nv)
      ctx%m=m; ctx%nv=nv; ctx%mu=mu; ctx%r=r
      allocate(ctx%diag(nv),ctx%block(size(prob%block))); ctx%diag=0.0_dp
      g=0.0_dp; sol%schur_assemblies=sol%schur_assemblies+1
      g(1:m)=prob%b; g(nv)=-penalty
      phi=dot_product(prob%b,z(1:m))-penalty*r
      if (r<=0.0_dp .or. .not.ieee_is_finite(r)) then; ok=.false.; return; end if
      phi=phi+mu*log(r); g(nv)=g(nv)+mu/r; ctx%diag(nv)=mu/(r*r)

      do k=1,size(prob%block)
         n=prob%block(k)%n; ctx%block(k)%category=prob%block(k)%category; ctx%block(k)%n=n
         select case(prob%block(k)%category)
         case(dsdp_sdp_block)
            allocate(s(n,n)); s=0.0_dp
            call data_add_scaled(prob%block(k),0,1.0_dp,s)
            do i=1,m; call data_add_scaled(prob%block(k),i,-z(i),s); end do
            do i=1,n; s(i,i)=s(i,i)+r; end do
            call spd_inverse_logdet(s,sinv,logdet,pd)
            if (.not.pd) then; ok=.false.; return; end if
            phi=phi+mu*logdet
            do i=1,m
               g(i)=g(i)-mu*data_dot(prob%block(k),i,sinv)
               v=data_schur_pair(prob%block(k),i,i,sinv)
               ctx%diag(i)=ctx%diag(i)+mu*v
               storage=data_storage(prob%block(k),i)
               if (storage==dsdp_data_dense) then
                  sol%dense_pair_evals=sol%dense_pair_evals+1
               else if (storage==dsdp_data_lowrank) then
                  sol%lowrank_pair_evals=sol%lowrank_pair_evals+1
               else
                  sol%sparse_pair_evals=sol%sparse_pair_evals+1
               end if
            end do
            g(nv)=g(nv)+mu*trace_mat(sinv)
            ctx%diag(nv)=ctx%diag(nv)+mu*sum(sinv*sinv)
            allocate(ctx%block(k)%sinv(n,n)); ctx%block(k)%sinv=sinv
            deallocate(s,sinv)

         case(dsdp_lp_block)
            allocate(sl(n),w(n))
            sl=prob%block(k)%cdiag-matmul(prob%block(k)%adiag,z(1:m))+r
            if (any(sl<=0.0_dp) .or. any(.not.ieee_is_finite(sl))) then; ok=.false.; return; end if
            phi=phi+mu*sum(log(sl)); w=1.0_dp/sl
            do i=1,m
               g(i)=g(i)-mu*dot_product(prob%block(k)%adiag(:,i),w)
            end do
            g(nv)=g(nv)+mu*sum(w); w=w*w
            do i=1,m
               ctx%diag(i)=ctx%diag(i)+mu*dot_product(prob%block(k)%adiag(:,i)*w, &
                  prob%block(k)%adiag(:,i))
            end do
            ctx%diag(nv)=ctx%diag(nv)+mu*sum(w)
            allocate(ctx%block(k)%w2(n)); ctx%block(k)%w2=w
            deallocate(sl,w)
         end select
      end do
      ctx%diag=max(ctx%diag,sqrt(tiny(1.0_dp)))
      ok=all(ieee_is_finite(g)) .and. all(ieee_is_finite(ctx%diag)) .and. ieee_is_finite(phi)
   end subroutine evaluate_newton_matrixfree

   subroutine schur_matvec(prob,ctx,p,q,reg,sol)
      type(dsdp_problem), intent(in) :: prob
      type(newton_context), intent(in) :: ctx
      real(dp), intent(in) :: p(:),reg
      real(dp), intent(out) :: q(:)
      type(dsdp_solution), intent(inout) :: sol
      real(dp), allocatable :: bmat(:,:),t(:,:),u(:)
      integer :: k,i,n,m,nv

      m=ctx%m; nv=ctx%nv; q=0.0_dp
      q(nv)=q(nv)+ctx%mu*p(nv)/(ctx%r*ctx%r)
      do k=1,size(prob%block)
         n=prob%block(k)%n
         select case(prob%block(k)%category)
         case(dsdp_sdp_block)
            allocate(bmat(n,n)); bmat=0.0_dp
            do i=1,m
               call data_add_scaled(prob%block(k),i,p(i),bmat)
            end do
            do i=1,n
               bmat(i,i)=bmat(i,i)-p(nv)
            end do
            t=matmul(matmul(ctx%block(k)%sinv,bmat),ctx%block(k)%sinv)
            do i=1,m
               q(i)=q(i)+ctx%mu*data_dot(prob%block(k),i,t)
            end do
            q(nv)=q(nv)-ctx%mu*trace_mat(t)
            deallocate(bmat,t)
         case(dsdp_lp_block)
            allocate(u(n))
            u=matmul(prob%block(k)%adiag,p(1:m))-p(nv)
            u=ctx%block(k)%w2*u
            q(1:m)=q(1:m)+ctx%mu*matmul(transpose(prob%block(k)%adiag),u)
            q(nv)=q(nv)-ctx%mu*sum(u)
            deallocate(u)
         end select
      end do
      if (reg>0.0_dp) q=q+reg*p
      sol%matrix_free_matvecs=sol%matrix_free_matvecs+1
   end subroutine schur_matvec

   subroutine pcg_matrixfree(prob,ctx,b,x,ctrl,sol,ok)
      type(dsdp_problem), intent(in) :: prob
      type(newton_context), intent(in) :: ctx
      real(dp), intent(in) :: b(:)
      real(dp), intent(out) :: x(:)
      type(dsdp_control), intent(in) :: ctrl
      type(dsdp_solution), intent(inout) :: sol
      logical, intent(out) :: ok
      real(dp), allocatable :: r(:),z(:),p(:),ap(:),diag(:)
      real(dp) :: rz,rzold,alpha,beta,pap,bnorm,rnorm,reg
      integer :: i,n,iters

      n=size(b); allocate(r(n),z(n),p(n),ap(n),diag(n)); x=0.0_dp
      reg=max(0.0_dp,ctrl%schur_regularization)
      diag=max(ctx%diag+reg,sqrt(tiny(1.0_dp)))
      r=b; z=r/diag; p=z; rz=dot_product(r,z); bnorm=norm2_dp(b)
      rnorm=norm2_dp(r); iters=0
      sol%cg_solves=sol%cg_solves+1; sol%matrix_free_cg_solves=sol%matrix_free_cg_solves+1
      if (rnorm<=ctrl%cg_tol*(1.0_dp+bnorm)) then; ok=.true.; return; end if
      do i=1,ctrl%cg_maxiter
         call schur_matvec(prob,ctx,p,ap,reg,sol)
         pap=dot_product(p,ap)
         if (pap<=tiny(1.0_dp) .or. .not.ieee_is_finite(pap)) exit
         alpha=rz/pap; x=x+alpha*p; r=r-alpha*ap; rnorm=norm2_dp(r); iters=i
         if (rnorm<=ctrl%cg_tol*(1.0_dp+bnorm)) then; ok=.true.; exit; end if
         z=r/diag; rzold=rz; rz=dot_product(r,z)
         if (rzold<=tiny(1.0_dp) .or. .not.ieee_is_finite(rz)) exit
         beta=rz/rzold; p=z+beta*p
      end do
      sol%cg_iterations=sol%cg_iterations+iters
      if (iters==0 .or. rnorm>ctrl%cg_tol*(1.0_dp+bnorm)) then
         ok=rnorm<=max(100.0_dp*ctrl%cg_tol*(1.0_dp+bnorm),1.0e-9_dp*(1.0_dp+bnorm))
      end if
   end subroutine pcg_matrixfree

   subroutine solve_schur_system(h,b,x,ctrl,sol,cache,ok)
      real(dp), intent(in) :: h(:,:),b(:)
      real(dp), intent(out) :: x(:)
      type(dsdp_control), intent(in) :: ctrl
      type(dsdp_solution), intent(inout) :: sol
      type(sparse_ldlt_cache), intent(inout) :: cache
      logical, intent(out) :: ok
      logical :: trycg,trysparse
      integer :: it,n
      real(dp) :: density

      n=size(b); density=schur_density(h,ctrl%sparse_schur_drop_tol)
      trysparse=ctrl%use_sparse_schur_factor .and. n>=max(1,ctrl%sparse_schur_threshold) .and. &
         density<=ctrl%sparse_schur_density_limit
      if (trysparse) then
         call sparse_ldlt_solve(h,b,x,ctrl%schur_regularization,ctrl%sparse_schur_drop_tol,cache,ok)
         sol%sparse_factor_solves=sol%sparse_factor_solves+1
         sol%sparse_symbolic_analyses=cache%symbolic_analyses
         sol%sparse_numeric_factorizations=cache%numeric_factorizations
         sol%schur_matrix_nnz=cache%matrix_nnz
         sol%schur_factor_nnz=cache%factor_nnz
         if (ok) return
         sol%sparse_factor_fallbacks=sol%sparse_factor_fallbacks+1
      end if

      trycg=ctrl%use_cg .and. .not.ctrl%cg_matrix_free .and. n>=max(1,ctrl%cg_threshold)
      if (trycg) then
         call pcg_spd(h,b,x,ctrl%schur_regularization,ctrl%cg_tol,ctrl%cg_maxiter,it,ok)
         sol%cg_solves=sol%cg_solves+1; sol%cg_iterations=sol%cg_iterations+it
         if (ok) return
         if (.not.ctrl%cg_fallback_direct) return
      end if
      call solve_spd(h,b,x,ctrl%schur_regularization,ok)
      sol%direct_schur_solves=sol%direct_schur_solves+1
   end subroutine solve_schur_system

   real(dp) function schur_density(h,drop_tol) result(density)
      real(dp), intent(in) :: h(:,:),drop_tol
      integer :: i,j,n,nnz
      real(dp) :: scale,tol
      n=size(h,1); nnz=0
      do j=1,n
         do i=1,n
            scale=sqrt(max(abs(h(i,i))*abs(h(j,j)),1.0_dp))
            tol=max(0.0_dp,drop_tol)*scale
            if (i==j .or. abs(h(i,j))>tol) nnz=nnz+1
         end do
      end do
      density=real(nnz,dp)/max(1.0_dp,real(n,dp)*real(n,dp))
   end function schur_density

   subroutine pcg_spd(a,b,x,reg,tol,maxiter,iters,ok)
      real(dp), intent(in) :: a(:,:),b(:),reg,tol
      real(dp), intent(out) :: x(:)
      integer, intent(in) :: maxiter
      integer, intent(out) :: iters
      logical, intent(out) :: ok
      real(dp), allocatable :: r(:),z(:),p(:),ap(:),diag(:)
      real(dp) :: rz,rzold,alpha,beta,pap,bnorm,rnorm
      integer :: i,n
      n=size(b); allocate(r(n),z(n),p(n),ap(n),diag(n)); x=0.0_dp
      do i=1,n; diag(i)=max(abs(a(i,i))+max(0.0_dp,reg),sqrt(tiny(1.0_dp))); end do
      r=b; z=r/diag; p=z; rz=dot_product(r,z); bnorm=norm2_dp(b)
      rnorm=norm2_dp(r); iters=0
      if (rnorm<=tol*(1.0_dp+bnorm)) then; ok=.true.; return; end if
      do i=1,maxiter
         ap=matmul(a,p)+max(0.0_dp,reg)*p
         pap=dot_product(p,ap)
         if (pap<=tiny(1.0_dp) .or. .not.ieee_is_finite(pap)) exit
         alpha=rz/pap; x=x+alpha*p; r=r-alpha*ap; rnorm=norm2_dp(r); iters=i
         if (rnorm<=tol*(1.0_dp+bnorm)) then; ok=.true.; return; end if
         z=r/diag; rzold=rz; rz=dot_product(r,z)
         if (rzold<=tiny(1.0_dp) .or. .not.ieee_is_finite(rz)) exit
         beta=rz/rzold; p=z+beta*p
      end do
      ok=rnorm<=max(100.0_dp*tol*(1.0_dp+bnorm),1.0e-9_dp*(1.0_dp+bnorm))
   end subroutine pcg_spd

   subroutine state_value(prob,z,mu,penalty,phi,ok)
      type(dsdp_problem), intent(in) :: prob
      real(dp), intent(in) :: z(:),mu,penalty
      real(dp), intent(out) :: phi
      logical, intent(out) :: ok
      real(dp), allocatable :: s(:,:),sl(:)
      real(dp) :: logdet,r
      integer :: k,i,n,m
      logical :: pd
      m=prob%m; r=z(m+1)
      if (r<=0.0_dp .or. .not.ieee_is_finite(r)) then; ok=.false.; phi=-huge(1.0_dp); return; end if
      phi=dot_product(prob%b,z(1:m))-penalty*r+mu*log(r)
      do k=1,size(prob%block)
         n=prob%block(k)%n
         if (prob%block(k)%category==dsdp_sdp_block) then
            allocate(s(n,n)); s=0.0_dp; call data_add_scaled(prob%block(k),0,1.0_dp,s)
            do i=1,m; call data_add_scaled(prob%block(k),i,-z(i),s); end do
            do i=1,n; s(i,i)=s(i,i)+r; end do
            call spd_logdet(s,logdet,pd)
            if (.not.pd) then; ok=.false.; phi=-huge(1.0_dp); return; end if
            phi=phi+mu*logdet; deallocate(s)
         else
            allocate(sl(n)); sl=prob%block(k)%cdiag-matmul(prob%block(k)%adiag,z(1:m))+r
            if (any(sl<=0.0_dp)) then; ok=.false.; phi=-huge(1.0_dp); return; end if
            phi=phi+mu*sum(log(sl)); deallocate(sl)
         end if
      end do
      ok=ieee_is_finite(phi)
   end subroutine state_value

   real(dp) function initial_shift(prob,y) result(r)
      type(dsdp_problem), intent(in) :: prob
      real(dp), intent(in) :: y(:)
      real(dp), allocatable :: s(:,:),sl(:)
      real(dp) :: dummy
      integer :: attempt,k,i,n
      logical :: ok,allok
      r=max(1.0e-8_dp,1.0e-8_dp*(1.0_dp+data_c_norm(prob)))
      do attempt=1,30
         allok=.true.
         do k=1,size(prob%block)
            n=prob%block(k)%n
            if (prob%block(k)%category==dsdp_sdp_block) then
               allocate(s(n,n)); s=0.0_dp; call data_add_scaled(prob%block(k),0,1.0_dp,s)
               do i=1,prob%m; call data_add_scaled(prob%block(k),i,-y(i),s); end do
               do i=1,n; s(i,i)=s(i,i)+r; end do
               call spd_logdet(s,dummy,ok); if (.not.ok) allok=.false.; deallocate(s)
            else
               allocate(sl(n)); sl=prob%block(k)%cdiag-matmul(prob%block(k)%adiag,y)+r
               if (any(sl<=0.0_dp)) allok=.false.; deallocate(sl)
            end if
         end do
         if (allok) return
         r=10.0_dp*r
      end do
      error stop 'dsdp_solve: unable to construct a strictly feasible shifted start'
   end function initial_shift

   subroutine recover_primal_and_metrics(prob,y,r,mu,sol,ok)
      type(dsdp_problem), intent(in) :: prob
      real(dp), intent(in) :: y(:),r,mu
      type(dsdp_solution), intent(inout) :: sol
      logical, intent(out) :: ok
      real(dp), allocatable :: s(:,:),sinv(:,:),sl(:),ax(:),rhs(:),hc(:,:),d(:),ad(:,:)
      real(dp) :: logdet,pobj,wmin,mineig,trv
      integer :: k,i,j,n,m
      logical :: pd,eigok,solveok
      m=prob%m
      if (.not.allocated(sol%x)) allocate(sol%x(size(prob%block)))
      allocate(ax(m),rhs(m),hc(m,m),d(m)); ax=0.0_dp; hc=0.0_dp

      do k=1,size(prob%block)
         n=prob%block(k)%n; sol%x(k)%category=prob%block(k)%category; sol%x(k)%n=n
         if (prob%block(k)%category==dsdp_sdp_block) then
            allocate(s(n,n)); s=0.0_dp; call data_add_scaled(prob%block(k),0,1.0_dp,s)
            do i=1,m; call data_add_scaled(prob%block(k),i,-y(i),s); end do
            do i=1,n; s(i,i)=s(i,i)+r; end do
            call spd_inverse_logdet(s,sinv,logdet,pd)
            if (.not.pd) then; ok=.false.; return; end if
            if (allocated(sol%x(k)%x)) deallocate(sol%x(k)%x)
            allocate(sol%x(k)%x(n,n)); sol%x(k)%x=mu*sinv
            do i=1,m; ax(i)=ax(i)+data_dot(prob%block(k),i,sol%x(k)%x); end do
            call add_sdp_hessian(prob%block(k),m,sinv,mu,hc,sol)
            deallocate(s,sinv)
         else
            allocate(sl(n)); sl=prob%block(k)%cdiag-matmul(prob%block(k)%adiag,y)+r
            if (any(sl<=0.0_dp)) then; ok=.false.; return; end if
            if (allocated(sol%x(k)%xdiag)) deallocate(sol%x(k)%xdiag)
            allocate(sol%x(k)%xdiag(n)); sol%x(k)%xdiag=mu/sl
            ax=ax+matmul(transpose(prob%block(k)%adiag),sol%x(k)%xdiag)
            do j=1,m; do i=1,j
               trv=mu*dot_product(prob%block(k)%adiag(:,i)/(sl*sl),prob%block(k)%adiag(:,j))
               hc(i,j)=hc(i,j)+trv; if (i/=j) hc(j,i)=hc(i,j)
            end do; end do
            deallocate(sl)
         end if
      end do

      rhs=prob%b-ax; call solve_spd(hc,rhs,d,1.0e-14_dp,solveok)
      if (.not.solveok) d=0.0_dp

      ax=0.0_dp; pobj=0.0_dp; mineig=huge(1.0_dp)
      do k=1,size(prob%block)
         n=prob%block(k)%n
         if (prob%block(k)%category==dsdp_sdp_block) then
            allocate(s(n,n)); s=0.0_dp; call data_add_scaled(prob%block(k),0,1.0_dp,s)
            do i=1,m; call data_add_scaled(prob%block(k),i,-y(i),s); end do
            call min_eigenvalue_sym(s,wmin,eigok); if (eigok) mineig=min(mineig,wmin)
            do i=1,n; s(i,i)=s(i,i)+r; end do
            call spd_inverse_logdet(s,sinv,logdet,pd)
            if (.not.pd) then; ok=.false.; return; end if
            allocate(ad(n,n)); ad=0.0_dp
            do i=1,m; call data_add_scaled(prob%block(k),i,d(i),ad); end do
            sol%x(k)%x=sol%x(k)%x+mu*matmul(matmul(sinv,ad),sinv)
            sol%x(k)%x=0.5_dp*(sol%x(k)%x+transpose(sol%x(k)%x))
            do i=1,m; ax(i)=ax(i)+data_dot(prob%block(k),i,sol%x(k)%x); end do
            pobj=pobj+data_dot(prob%block(k),0,sol%x(k)%x)
            deallocate(s,sinv,ad)
         else
            allocate(sl(n)); sl=prob%block(k)%cdiag-matmul(prob%block(k)%adiag,y)
            mineig=min(mineig,minval(sl)); sl=sl+r
            sol%x(k)%xdiag=sol%x(k)%xdiag+mu*matmul(prob%block(k)%adiag,d)/(sl*sl)
            where (sol%x(k)%xdiag<0.0_dp .and. sol%x(k)%xdiag>-1.0e-12_dp) sol%x(k)%xdiag=0.0_dp
            ax=ax+matmul(transpose(prob%block(k)%adiag),sol%x(k)%xdiag)
            pobj=pobj+dot_product(prob%block(k)%cdiag,sol%x(k)%xdiag); deallocate(sl)
         end if
      end do
      sol%pobj=pobj; sol%dobj=dot_product(prob%b,y)
      sol%pinfeas=norm2_dp(ax-prob%b)/(1.0_dp+norm2_dp(prob%b))
      sol%relgap=abs(sol%pobj-sol%dobj)/(1.0_dp+abs(sol%pobj)+abs(sol%dobj))
      sol%min_dual_eig=mineig
      ok=ieee_is_finite(sol%pobj) .and. ieee_is_finite(sol%dobj) .and. &
         ieee_is_finite(sol%pinfeas) .and. ieee_is_finite(sol%relgap)
   end subroutine recover_primal_and_metrics

   pure real(dp) function trace_mat(a) result(v)
      real(dp), intent(in) :: a(:,:)
      integer :: i
      v=0.0_dp; do i=1,min(size(a,1),size(a,2)); v=v+a(i,i); end do
   end function trace_mat

   real(dp) function data_c_norm(prob) result(v)
      type(dsdp_problem), intent(in) :: prob
      integer :: k
      v=0.0_dp
      do k=1,size(prob%block)
         if (prob%block(k)%category==dsdp_sdp_block) then; v=v+data_fnorm2(prob%block(k),0)
         else; v=v+dot_product(prob%block(k)%cdiag,prob%block(k)%cdiag); end if
      end do
      v=sqrt(max(0.0_dp,v))
   end function data_c_norm

   integer function count_problem_nnz(prob) result(nnz)
      type(dsdp_problem), intent(in) :: prob
      integer :: k,i
      nnz=0
      do k=1,size(prob%block)
         if (prob%block(k)%category==dsdp_sdp_block) then
            nnz=nnz+data_nnz(prob%block(k),0)
            do i=1,prob%m; nnz=nnz+data_nnz(prob%block(k),i); end do
         else
            nnz=nnz+count(abs(prob%block(k)%cdiag)>tiny(1.0_dp))
            nnz=nnz+count(abs(prob%block(k)%adiag)>tiny(1.0_dp))
         end if
      end do
   end function count_problem_nnz

end module rdsdp_solver
