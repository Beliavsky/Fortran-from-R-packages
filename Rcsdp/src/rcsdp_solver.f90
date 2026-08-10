! Predictor-corrector SDP solver translated from CSDP 6.1.1.
! The Newton equations and stopping measures follow CSDP's sdp.c/easysdp.c.
! v0.2.0 added sparse Schur kernels; v0.3.0 adds fill products, scaling, and Lanczos line search.
! See LICENSE (CPL-1.0).
module rcsdp_solver
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rcsdp_kinds, only : dp
   use rcsdp_types, only : csdp_problem, csdp_control, csdp_solution, csdp_block_matrix, &
      csdp_success, csdp_primal_infeas, csdp_dual_infeas, csdp_partial_success, &
      csdp_max_iterations, csdp_primal_edge, csdp_dual_edge, csdp_no_progress, &
      csdp_singular, csdp_numeric_failure
   use rcsdp_block_ops, only : zero_mat, make_i, add_mat, addscaledmat, sym_mat, mat_mult, &
      trace_prod, fnorm, chol, inverse_from_chol, line_search_pd, calc_pobj, block_is_pd
   use rcsdp_problem_mod, only : build_dense_constraints, dense_constraint, op_a, op_at, op_o, validate_problem
   use rcsdp_sparse_ops, only : csdp_sparse_workspace, build_sparse_workspace, op_a_sparse, op_at_sparse, &
      op_o_sparse, apply_o_sparse
   use rcsdp_fill_ops, only : csdp_fill_workspace, build_fill_workspace, mat_multspa, mat_multspb, &
      mat_multspc, apply_o_fill
   use rcsdp_linalg, only : norm2_dp, potrf_upper, solve_spd_factored
   implicit none
   private
   public :: csdp, initsoln, pinfeas, dinfeas, calc_dobj

contains

   pure real(dp) function calc_dobj(prob,y) result(v)
      type(csdp_problem), intent(in) :: prob
      real(dp), intent(in) :: y(:)
      v = dot_product(prob%b,y) + prob%constant_offset
   end function calc_dobj

   real(dp) function pinfeas(prob,x,adense) result(v)
      type(csdp_problem), intent(in) :: prob
      type(csdp_block_matrix), intent(in) :: x
      type(csdp_block_matrix), intent(in), optional :: adense(:)
      real(dp), allocatable :: r(:)
      allocate(r(size(prob%b)))
      if (present(adense)) then
         call op_a(prob,x,r,adense)
      else
         call op_a_sparse(prob,x,r)
      end if
      r=r-prob%b
      v=norm2_dp(r)/(1.0_dp+norm2_dp(prob%b))
   end function pinfeas

   real(dp) function dinfeas(prob,y,z,adense) result(v)
      type(csdp_problem), intent(in) :: prob
      real(dp), intent(in) :: y(:)
      type(csdp_block_matrix), intent(in) :: z
      type(csdp_block_matrix), intent(in), optional :: adense(:)
      type(csdp_block_matrix) :: w
      if (present(adense)) then
         call op_at(prob,y,w,adense)
      else
         call op_at_sparse(prob,y,w)
      end if
      call addscaledmat(w,-1.0_dp,prob%c,w)
      call addscaledmat(w,-1.0_dp,z,w)
      v=fnorm(w)/(1.0_dp+fnorm(prob%c))
   end function dinfeas

   subroutine initsoln(prob,x,y,z,adense)
      type(csdp_problem), intent(in) :: prob
      type(csdp_block_matrix), intent(out) :: x,z
      real(dp), allocatable, intent(out) :: y(:)
      type(csdp_block_matrix), intent(in), optional :: adense(:)
      type(csdp_block_matrix) :: ai
      real(dp) :: alpha,beta,maxnrma,nrma,nrmc
      integer :: i,k,n
      n=prob%order()
      alpha=0.0_dp
      maxnrma=0.0_dp
      do i=1,size(prob%b)
         if (present(adense)) then
            nrma=fnorm(adense(i))
         else
            call dense_constraint(prob,i,ai)
            nrma=fnorm(ai)
         end if
         maxnrma=max(maxnrma,nrma)
         alpha=max(alpha,(1.0_dp+abs(prob%b(i)))/(1.0_dp+nrma))
      end do
      alpha=max(1.0_dp,real(n,dp)*alpha)
      nrmc=fnorm(prob%c)
      beta=(1.0_dp+max(nrmc,maxnrma))/sqrt(max(1.0_dp,real(n,dp)))
      x=prob%c
      z=prob%c
      call make_i(x)
      call make_i(z)
      do k=1,size(x%block)
         if (allocated(x%block(k)%diag)) then
            x%block(k)%diag=10.0_dp*alpha*x%block(k)%diag
            z%block(k)%diag=10.0_dp*beta*z%block(k)%diag
         else
            x%block(k)%mat=10.0_dp*alpha*x%block(k)%mat
            z%block(k)%mat=10.0_dp*beta*z%block(k)%mat
         end if
      end do
      allocate(y(size(prob%b)))
      y=0.0_dp
   end subroutine initsoln

   subroutine csdp(prob,sol,control,x0,y0,z0)
      type(csdp_problem), intent(in) :: prob
      type(csdp_solution), intent(out) :: sol
      type(csdp_control), intent(in), optional :: control
      type(csdp_block_matrix), intent(in), optional :: x0,z0
      real(dp), intent(in), optional :: y0(:)

      type(csdp_control) :: ctrl
      type(csdp_block_matrix), allocatable :: adense(:)
      type(csdp_sparse_workspace) :: sparse_work
      type(csdp_fill_workspace) :: fill_work
      type(csdp_block_matrix) :: x,z,zi,zchol,rinv
      type(csdp_block_matrix) :: fd,aty,w1,w2,w3,dx,dz,dz1,dx1,ident
      type(csdp_block_matrix) :: bestx,bestz,trialx,trialz
      real(dp), allocatable :: y(:),besty(:),dy(:),dy1(:),rhs(:),fp(:),tmpv(:),o(:,:),ofac(:,:),schur_scale(:)
      real(dp) :: pobj,dobj,gap,relgap,relp,reld
      real(dp) :: bestmeasure,measure
      real(dp) :: alphap,alphad,alphap1,alphad1,minalpha,mystepfrac
      real(dp) :: muk,muplus,gamma,mu,diagnrm,diagadd
      real(dp) :: perturbfac,normc,norma,pinfmeas,dinfmeas,den
      real(dp) :: schur_seconds,tic,toc,max_diagadd
      integer :: iter,info,tries,n,m,status,stagnation
      integer :: spairs,dprods,total_spairs,total_dprods,schur_assemblies,total_refinements
      integer :: fill_sparse_products,fill_dense_products,lanczos_linesearches
      logical :: ok,use_fill
      character(len=:), allocatable :: msg

      ctrl=csdp_control()
      if (present(control)) ctrl=control
      call validate_problem(prob,ok,msg)
      if (.not.ok) error stop 'csdp: '//msg
      call build_sparse_workspace(prob,sparse_work)
      use_fill = ctrl%use_sparse_schur .and. ctrl%use_fill_products
      if (use_fill) call build_fill_workspace(prob,fill_work)
      if (.not.ctrl%use_sparse_schur) call build_dense_constraints(prob,adense)
      schur_seconds=0.0_dp
      total_spairs=0
      total_dprods=0
      schur_assemblies=0
      total_refinements=0
      fill_sparse_products=0
      fill_dense_products=0
      lanczos_linesearches=0
      max_diagadd=0.0_dp
      n=prob%order()
      m=prob%nconstraints()

      if (present(x0) .and. present(y0) .and. present(z0)) then
         x=x0; z=z0
         allocate(y(size(y0))); y=y0
      else
         if (ctrl%use_sparse_schur) then
            call initsoln(prob,x,y,z)
         else
            call initsoln(prob,x,y,z,adense)
         end if
      end if
      if (size(y)/=m) error stop 'csdp: initial y has wrong size'

      fd=prob%c; aty=prob%c; w1=prob%c; w2=prob%c; w3=prob%c
      dx=prob%c; dz=prob%c; dz1=prob%c; dx1=prob%c; ident=prob%c
      trialx=prob%c; trialz=prob%c
      call make_i(ident)
      allocate(dy(m),dy1(m),rhs(m),fp(m),tmpv(m),o(m,m),ofac(m,m),besty(m),schur_scale(m))

      normc=fnorm(prob%c)
      norma=norm2_dp(prob%b)
      perturbfac=0.0_dp
      if (ctrl%perturbobj .and. n>0) perturbfac=1.0e-6_dp*normc/sqrt(real(n,dp))

      ! Initial inverse Z^{-1} and positive-definiteness checks.
      zchol=z
      info=chol(zchol)
      if (info/=0) then
         call finish(csdp_singular,0)
         return
      end if
      call inverse_from_chol(zchol,zi,rinv,info)
      if (info/=0) then
         call finish(csdp_singular,0)
         return
      end if
      w1=x
      info=chol(w1)
      if (info/=0) then
         call finish(csdp_singular,0)
         return
      end if

      bestx=x; bestz=z; besty=y
      bestmeasure=huge(1.0_dp)
      alphap=0.0_dp; alphad=0.0_dp
      mu=huge(1.0_dp)
      status=csdp_max_iterations
      stagnation=0

      do iter=0,ctrl%maxiter
         call metrics(pobj,dobj,gap,relgap,relp,reld)
         measure=max(relgap/max(ctrl%objtol,tiny(1.0_dp)), &
                     relp/max(ctrl%axtol,tiny(1.0_dp)), &
                     reld/max(ctrl%atytol,tiny(1.0_dp)))
         if (measure<bestmeasure .and. all_finite_metrics(pobj,dobj,relgap,relp,reld)) then
            bestmeasure=measure; bestx=x; bestz=z; besty=y; stagnation=0
         else
            stagnation=stagnation+1
         end if

         if (ctrl%printlevel>=1) then
            write(*,'("Iter: ",i3," Ap: ",es9.2," Pobj: ",es15.7," Ad: ",es9.2," Dobj: ",es15.7)') &
               iter,alphap,pobj,alphad,dobj
         end if

         if (relgap<=ctrl%objtol .and. relp<=ctrl%axtol .and. reld<=ctrl%atytol) then
            status=csdp_success
            exit
         end if
         if (.not.all_finite_metrics(pobj,dobj,relgap,relp,reld)) then
            status=csdp_numeric_failure
            exit
         end if
         if (iter==ctrl%maxiter) then
            status=csdp_max_iterations
            exit
         end if
         if (stagnation>60 .and. bestmeasure>1.0_dp) then
            status=csdp_no_progress
            exit
         end if

         ! CSDP infeasibility certificates.
         call apply_at_local(y,w1)
         call addscaledmat(w1,-1.0_dp,z,w1)
         den=fnorm(w1)
         if (den>0.0_dp) then
            pinfmeas=-dobj/den
            if (pinfmeas>ctrl%pinftol .and. relp>ctrl%axtol) then
               if (abs(dobj)>tiny(1.0_dp)) then
                  y=-y/dobj
                  call zero_mat(w1)
                  call addscaledmat(w1,-1.0_dp/dobj,z,z)
               end if
               status=csdp_primal_infeas
               bestx=x; bestz=z; besty=y
               exit
            end if
         end if
         call apply_a_local(x,tmpv)
         den=norm2_dp(tmpv)
         if (den>0.0_dp) then
            dinfmeas=trace_prod(prob%c,x)/den
            if (dinfmeas>ctrl%dinftol .and. reld>ctrl%atytol) then
               den=trace_prod(prob%c,x)
               if (abs(den)>tiny(1.0_dp)) then
                  call zero_mat(w1)
                  call addscaledmat(w1,1.0_dp/den,x,x)
               end if
               status=csdp_dual_infeas
               bestx=x; bestz=z; besty=y
               exit
            end if
         end if

         if (iter>1) then
            minalpha=min(alphap,alphad)
         else
            minalpha=1.0_dp
         end if
         mystepfrac=ctrl%minstepfrac+minalpha*(ctrl%maxstepfrac-ctrl%minstepfrac)

         ! O(.) = A(Z^{-1} A'(. ) X).
         call assemble_schur_local(zi,x,o)
         call factor_schur(o,ofac,schur_scale,info,diagnrm,diagadd)
         max_diagadd=max(max_diagadd,diagadd)
         if (info/=0) then
            status=csdp_singular
            exit
         end if

         ! Fd = Z + C - A'(y), with the same small objective perturbation idea
         ! used by CSDP while far from convergence.
         call apply_at_local(y,aty)
         call addscaledmat(z,1.0_dp,prob%c,fd)
         if (perturbfac>0.0_dp .and. bestmeasure>1.0e3_dp) then
            call addscaledmat(fd,-perturbfac,ident,fd)
         else if (perturbfac>0.0_dp .and. bestmeasure<1.0e3_dp) then
            call addscaledmat(fd,-perturbfac*(max(bestmeasure,0.0_dp)/1000.0_dp)**1.5_dp,ident,fd)
         end if
         call addscaledmat(fd,-1.0_dp,aty,fd)

         ! Affine predictor rhs = -b + A(Z^{-1} Fd X).
         if (use_fill) then
            call mat_multspb(1.0_dp,0.0_dp,zi,fd,w1,fill_work,ctrl%fill_density_limit, &
               fill_sparse_products,fill_dense_products)
            call mat_multspc(1.0_dp,0.0_dp,w1,x,w2,fill_work,ctrl%fill_density_limit, &
               fill_sparse_products,fill_dense_products)
         else
            call mat_mult(1.0_dp,0.0_dp,zi,fd,w1)
            call mat_mult(1.0_dp,0.0_dp,w1,x,w2)
         end if
         call apply_a_local(w2,rhs)
         rhs=rhs-prob%b
         dy=rhs*schur_scale
         call solve_spd_factored(ofac,dy,info)
         dy=dy*schur_scale
         if (info/=0) then
            status=csdp_singular
            exit
         end if
         if (.not.ctrl%fastmode .and. iter>1) call refine_schur_local(rhs,dy,ofac)

         ! Predictor directions.  With fill products enabled, preserve the
         ! already-computed Zi*Fd from the affine RHS and mirror CSDP's
         ! original construction dX = -(I - Zi*Fd + Zi*A'(dy))*X.
         if (use_fill) then
            call apply_at_local(dy,w3)
            dz=w3
            call addscaledmat(dz,-1.0_dp,fd,dz)
            call mat_multspb(1.0_dp,0.0_dp,zi,w3,w2,fill_work,ctrl%fill_density_limit, &
               fill_sparse_products,fill_dense_products)
            w3=ident
            call addscaledmat(w3,-1.0_dp,w1,w3)
            call addscaledmat(w3,1.0_dp,w2,w3)
            call mat_mult(-1.0_dp,0.0_dp,w3,x,dx)
         else
            call apply_at_local(dy,dz)
            call addscaledmat(dz,-1.0_dp,fd,dz)
            call apply_at_local(dy,w3)
            call mat_mult(1.0_dp,0.0_dp,zi,fd,w1)
            call mat_mult(1.0_dp,0.0_dp,w1,x,w2)
            call mat_mult(1.0_dp,0.0_dp,zi,w3,w1)
            call mat_mult(1.0_dp,0.0_dp,w1,x,dx)
            call addscaledmat(w2,-1.0_dp,dx,dx)
            call addscaledmat(dx,-1.0_dp,x,dx)
         end if
         call sym_mat(dx)

         alphap1=line_search_local(dx,x,mystepfrac,1.0_dp)
         alphad1=line_search_local(dz,z,mystepfrac,1.0_dp)
         if (alphap1<0.0_dp .or. alphad1<0.0_dp) then
            status=csdp_numeric_failure
            exit
         end if
         call addscaledmat(x,alphap1,dx,trialx)
         call addscaledmat(z,alphad1,dz,trialz)

         muk=abs(trace_prod(x,z))/max(1.0_dp,real(n,dp))
         muplus=trace_prod(trialx,trialz)/max(1.0_dp,real(n,dp))
         if (muplus<0.0_dp) muplus=0.5_dp*muk
         if (muk<=tiny(1.0_dp)) then
            mu=0.0_dp
         else
            gamma=max(0.0_dp,muplus/muk)
            if (relp<0.1_dp*ctrl%axtol .and. reld<0.1_dp*ctrl%atytol .and. &
                alphap>0.2_dp .and. alphad>0.2_dp .and. mu>1.0e-6_dp .and. &
                gamma<1.0_dp .and. alphap+alphad>1.0_dp) then
               mu=muk*gamma**(alphap+alphad)
            else
               mu=min(muplus,0.9_dp*muk)
            end if
            if (relp<0.9_dp*ctrl%axtol .and. reld<0.9_dp*ctrl%atytol) mu=min(mu,0.5_dp*muk)
         end if
         if (ctrl%affine) mu=0.0_dp

         ! Corrector: Fp = b - A(X+dX).
         call addscaledmat(x,1.0_dp,dx,w1)
         call apply_a_local(w1,fp)
         fp=prob%b-fp

         ! rhs = A[Z^{-1}(mu I - dZ dX)] - Fp.
         if (use_fill) then
            call mat_multspa(-1.0_dp,0.0_dp,dz,dx,w1,fill_work,ctrl%fill_density_limit, &
               fill_sparse_products,fill_dense_products)
            call addscaledmat(w1,mu,ident,w1)
            call mat_multspc(1.0_dp,0.0_dp,zi,w1,w2,fill_work,ctrl%fill_density_limit, &
               fill_sparse_products,fill_dense_products)
         else
            call mat_mult(-1.0_dp,0.0_dp,dz,dx,w1)
            call addscaledmat(w1,mu,ident,w1)
            call mat_mult(1.0_dp,0.0_dp,zi,w1,w2)
         end if
         call apply_a_local(w2,rhs)
         rhs=rhs-fp
         dy1=rhs*schur_scale
         call solve_spd_factored(ofac,dy1,info)
         dy1=dy1*schur_scale
         if (info/=0) then
            status=csdp_singular
            exit
         end if
         if (.not.ctrl%fastmode .and. iter>1) call refine_schur_local(rhs,dy1,ofac)

         call apply_at_local(dy1,dz1)
         ! dX1 = Zi * (-dZ1*X - dZ*dX + mu*I).
         if (use_fill) then
            call mat_multspa(-1.0_dp,0.0_dp,dz,dx,w1,fill_work,ctrl%fill_density_limit, &
               fill_sparse_products,fill_dense_products)
            call addscaledmat(w1,mu,ident,w1)
            call mat_multspa(-1.0_dp,1.0_dp,dz1,x,w1,fill_work,ctrl%fill_density_limit, &
               fill_sparse_products,fill_dense_products)
            call mat_mult(1.0_dp,0.0_dp,zi,w1,dx1)
         else
            call mat_mult(-1.0_dp,0.0_dp,dz,dx,w1)
            call addscaledmat(w1,mu,ident,w1)
            call mat_mult(-1.0_dp,1.0_dp,dz1,x,w1)
            call mat_mult(1.0_dp,0.0_dp,zi,w1,dx1)
         end if
         call sym_mat(dx1)

         call add_mat(dx1,dx)
         call add_mat(dz1,dz)
         dy=dy+dy1

         alphap=line_search_local(dx,x,mystepfrac,1.0_dp)
         alphad=line_search_local(dz,z,mystepfrac,1.0_dp)
         if (alphap<0.0_dp .or. alphad<0.0_dp) then
            status=csdp_numeric_failure
            exit
         end if
         if (alphap<ctrl%minstepp) then
            status=csdp_primal_edge
            exit
         end if
         if (alphad<ctrl%minstepd) then
            status=csdp_dual_edge
            exit
         end if

         call addscaledmat(x,alphap,dx,trialx)
         tries=0
         do while (.not.block_is_pd(trialx) .and. tries<30)
            alphap=0.9_dp*alphap
            call addscaledmat(x,alphap,dx,trialx)
            tries=tries+1
         end do
         if (tries>=30) then
            status=csdp_primal_edge
            exit
         end if
         call addscaledmat(z,alphad,dz,trialz)
         tries=0
         do while (.not.block_is_pd(trialz) .and. tries<30)
            alphad=0.9_dp*alphad
            call addscaledmat(z,alphad,dz,trialz)
            tries=tries+1
         end do
         if (tries>=30) then
            status=csdp_dual_edge
            exit
         end if

         x=trialx; z=trialz; y=y+alphad*dy

         zchol=z
         info=chol(zchol)
         if (info/=0) then
            status=csdp_singular
            exit
         end if
         call inverse_from_chol(zchol,zi,rinv,info)
         if (info/=0) then
            status=csdp_singular
            exit
         end if
      end do

      if (status/=csdp_primal_infeas .and. status/=csdp_dual_infeas) then
         x=bestx; z=bestz; y=besty
         call metrics(pobj,dobj,gap,relgap,relp,reld)

         ! CSDP's optional negative objective-gap tweak: move y in the b
         ! direction and Z by A'(b), limited by positive definiteness of Z.
         if (ctrl%tweakgap .and. .not.ctrl%usexzgap .and. dobj-pobj<0.0_dp .and. norma>tiny(1.0_dp)) then
            dy=prob%b
            call apply_at_local(dy,dz)
            alphad=line_search_local(dz,z,1.0_dp,-(dobj-pobj)/(norma*norma))
            if (alphad>0.0_dp) then
               y=y+alphad*dy
               call addscaledmat(z,alphad,dz,z)
               call metrics(pobj,dobj,gap,relgap,relp,reld)
            end if
         end if

         measure=max(relgap/max(ctrl%objtol,tiny(1.0_dp)), &
                     relp/max(ctrl%axtol,tiny(1.0_dp)), &
                     reld/max(ctrl%atytol,tiny(1.0_dp)))
         if (measure<=1.0_dp) then
            status=csdp_success
         else if (measure<1000.0_dp .and. status/=csdp_success) then
            status=csdp_partial_success
         end if
      end if
      call finish(status,min(iter,ctrl%maxiter))

   contains

      subroutine apply_a_local(mat,vec)
         type(csdp_block_matrix), intent(in) :: mat
         real(dp), intent(out) :: vec(:)
         if (ctrl%use_sparse_schur) then
            call op_a_sparse(prob,mat,vec)
         else
            call op_a(prob,mat,vec,adense)
         end if
      end subroutine apply_a_local

      subroutine apply_at_local(vec,mat)
         real(dp), intent(in) :: vec(:)
         type(csdp_block_matrix), intent(out) :: mat
         if (ctrl%use_sparse_schur) then
            call op_at_sparse(prob,vec,mat)
         else
            call op_at(prob,vec,mat,adense)
         end if
      end subroutine apply_at_local

      subroutine assemble_schur_local(zinv,xmat,omat)
         type(csdp_block_matrix), intent(in) :: zinv,xmat
         real(dp), intent(out) :: omat(:,:)
         call cpu_time(tic)
         if (ctrl%use_sparse_schur) then
            call op_o_sparse(prob,zinv,xmat,omat,sparse_work,spairs,dprods)
            total_spairs=total_spairs+spairs
            total_dprods=total_dprods+dprods
         else
            call op_o(prob,zinv,xmat,omat,adense)
         end if
         call cpu_time(toc)
         schur_seconds=schur_seconds+max(0.0_dp,toc-tic)
         schur_assemblies=schur_assemblies+1
      end subroutine assemble_schur_local

      subroutine refine_schur_local(bvec,xvec,rchol)
         real(dp), intent(in) :: bvec(:),rchol(:,:)
         real(dp), intent(inout) :: xvec(:)
         real(dp), allocatable :: ax(:),resid(:),corr(:),trial(:),trialax(:)
         real(dp) :: besterr,err,denom
         integer :: kref,ifail,noimprove

         allocate(ax(size(bvec)),resid(size(bvec)),corr(size(bvec)), &
            trial(size(bvec)),trialax(size(bvec)))
         if (use_fill) then
            call apply_o_fill(prob,zi,x,xvec,ax,fill_work,ctrl%fill_density_limit, &
               fill_sparse_products,fill_dense_products)
         else
            call apply_o_sparse(prob,zi,x,xvec,ax)
         end if
         resid=bvec-ax
         denom=1.0_dp+norm2_dp(bvec)
         besterr=norm2_dp(resid)/denom
         if (besterr<=1.0e-14_dp) return
         noimprove=0
         do kref=1,20
            corr=resid*schur_scale
            call solve_spd_factored(rchol,corr,ifail)
            corr=corr*schur_scale
            if (ifail/=0) exit
            trial=xvec+corr
            if (use_fill) then
               call apply_o_fill(prob,zi,x,trial,trialax,fill_work,ctrl%fill_density_limit, &
                  fill_sparse_products,fill_dense_products)
            else
               call apply_o_sparse(prob,zi,x,trial,trialax)
            end if
            err=norm2_dp(bvec-trialax)/denom
            if (err<besterr) then
               xvec=trial
               ax=trialax
               resid=bvec-ax
               besterr=err
               total_refinements=total_refinements+1
               noimprove=0
            else
               noimprove=noimprove+1
            end if
            if (besterr<=1.0e-14_dp .or. noimprove>=3) exit
         end do
      end subroutine refine_schur_local

      subroutine metrics(pp,dd,gg,rg,rp,rd)
         real(dp), intent(out) :: pp,dd,gg,rg,rp,rd
         pp=calc_pobj(prob%c,x,prob%constant_offset)
         dd=calc_dobj(prob,y)
         if (ctrl%usexzgap) then
            gg=max(0.0_dp,trace_prod(x,z))
         else
            gg=max(0.0_dp,dd-pp)
         end if
         rg=gg/(1.0_dp+abs(pp)+abs(dd))
         if (ctrl%use_sparse_schur) then
            rp=pinfeas(prob,x)
            rd=dinfeas(prob,y,z)
         else
            rp=pinfeas(prob,x,adense)
            rd=dinfeas(prob,y,z,adense)
         end if
      end subroutine metrics

      logical function all_finite_metrics(pp,dd,rg,rp,rd) result(f)
         real(dp), intent(in) :: pp,dd,rg,rp,rd
         f=ieee_is_finite(pp) .and. ieee_is_finite(dd) .and. ieee_is_finite(rg) .and. &
            ieee_is_finite(rp) .and. ieee_is_finite(rd)
      end function all_finite_metrics

      subroutine factor_schur(a,r,scale,ifail,dn,da)
         real(dp), intent(in) :: a(:,:)
         real(dp), intent(out) :: r(:,:),scale(:)
         integer, intent(out) :: ifail
         real(dp), intent(out) :: dn,da
         real(dp), allocatable :: basea(:,:), diagv(:)
         real(dp) :: diagfact,unit,mindiag
         integer :: retries,jj,nn

         nn=size(a,1)
         allocate(basea(nn,nn),diagv(nn))
         basea=0.5_dp*(a+transpose(a))
         do jj=1,nn
            diagv(jj)=basea(jj,jj)
         end do
         dn=sqrt(dot_product(diagv,diagv))
         mindiag=minval(diagv)
         unit=1.0e-17_dp*dn/sqrt(max(1.0_dp,real(nn,dp)))
         if (unit<=tiny(1.0_dp)) unit=tiny(1.0_dp)
         retries=0
         diagfact=0.0_dp

         do
            da=unit*diagfact
            do while (mindiag+da<=0.0_dp)
               retries=retries+1
               if (abs(diagfact)<=tiny(1.0_dp)) then
                  diagfact=0.1_dp
               else
                  diagfact=10.0_dp*diagfact
               end if
               da=unit*diagfact
               if (retries>15) exit
            end do
            if (retries>15) then
               ifail=1
               return
            end if

            r=basea
            do jj=1,nn
               r(jj,jj)=r(jj,jj)+da
            end do
            if (ctrl%use_schur_scaling) then
               do jj=1,nn
                  scale(jj)=min(1.0e30_dp,1.0_dp/sqrt(max(r(jj,jj),tiny(1.0_dp))))
               end do
               r=r*spread(scale,2,nn)*spread(scale,1,nn)
            else
               scale=1.0_dp
            end if

            call potrf_upper(r,ifail)
            if (ifail==0) return
            if (retries>=15) return
            if (retries==0) diagfact=0.1_dp
            retries=retries+1
            diagfact=10.0_dp*diagfact
         end do
      end subroutine factor_schur

      real(dp) function line_search_local(dmat,xmat,stepfrac,startfrac) result(aout)
         type(csdp_block_matrix), intent(in) :: dmat,xmat
         real(dp), intent(in) :: stepfrac,startfrac
         logical :: used
         aout=line_search_pd(dmat,xmat,stepfrac,startfrac,ctrl%use_lanczos_linesearch, &
            ctrl%lanczos_threshold,ctrl%lanczos_iterations,used)
         if (used) lanczos_linesearches=lanczos_linesearches+1
      end function line_search_local

      subroutine finish(code,iters)
         integer, intent(in) :: code,iters
         real(dp) :: pp,dd,gg,rg,rp,rd
         if (.not.allocated(sol%y)) allocate(sol%y(size(y)))
         sol%x=x; sol%z=z; sol%y=y
         call metrics(pp,dd,gg,rg,rp,rd)
         sol%pobj=pp; sol%dobj=dd; sol%relgap=rg; sol%pinfeas=rp; sol%dinfeas=rd
         sol%status=code; sol%iterations=iters
         sol%constraint_nnz=sparse_work%constraint_nnz
         sol%sparse_constraint_blocks=sparse_work%sparse_blocks
         sol%dense_constraint_blocks=sparse_work%dense_blocks
         sol%schur_assemblies=schur_assemblies
         sol%schur_sparse_pairs=total_spairs
         sol%schur_dense_products=total_dprods
         sol%schur_refinements=total_refinements
         sol%schur_seconds=schur_seconds
         if (allocated(fill_work%block)) then
            sol%fill_nnz=fill_work%fill_nnz
            sol%fill_full_entries=fill_work%full_entries
         else
            sol%fill_nnz=0
            sol%fill_full_entries=0
         end if
         sol%fill_sparse_products=fill_sparse_products
         sol%fill_dense_products=fill_dense_products
         sol%lanczos_linesearches=lanczos_linesearches
         sol%schur_diagadd=max_diagadd
      end subroutine finish

   end subroutine csdp

end module rcsdp_solver
