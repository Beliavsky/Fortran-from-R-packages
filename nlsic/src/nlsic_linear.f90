! SPDX-License-Identifier: GPL-2.0-only
module nlsic_linear
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use nlsic_kinds, only : dp
   use nlsic_types, only : lsi_result, LSI_SUCCESS, LSI_INVALID_INPUT, &
      LSI_INFEASIBLE, LSI_RANK_DEFICIENT, LSI_NUMERICAL
   use nlsic_linalg, only : pseudoinverse, matrix_rank, null_space, &
      least_norm_ls, symmetric_eigen, dense_solve, vecnorm, identity_matrix
   use nlsic_nnls, only : ldp_solve
   implicit none
   private
   public :: ldp, ls_ln, ls_ln_svd, lsi, lsi_ln, lsie_ln, lsi_reg
   public :: tls, nulla, pnull, uplo_to_uco, vecnorm, ls_ln_multi
   interface ls_ln
      module procedure ls_ln_vector
   end interface ls_ln
   interface nulla
      module procedure nulla_matrix, nulla_vector
   end interface nulla
contains

   subroutine init_result(res,n,m)
      type(lsi_result), intent(out) :: res
      integer, intent(in) :: n,m
      allocate(res%x(max(0,n)),res%residuals(max(0,m)))
      res%x=0.0_dp; res%residuals=0.0_dp; res%rnorm=0.0_dp
      res%objective=0.0_dp; res%rank=0; res%status=LSI_SUCCESS; res%lambda=0.0_dp
   end subroutine init_result

   subroutine finish_ls(a,b,res)
      real(dp), intent(in) :: a(:,:),b(:)
      type(lsi_result), intent(inout) :: res
      res%residuals=matmul(a,res%x)-b
      res%rnorm=vecnorm(res%residuals)
      res%objective=res%rnorm*res%rnorm
   end subroutine finish_ls

   subroutine ldp(u,co,res,rcond)
      real(dp), intent(in) :: u(:,:),co(:)
      type(lsi_result), intent(out) :: res
      real(dp), intent(in), optional :: rcond
      logical :: feasible
      integer :: n
      n=size(u,2); call init_result(res,n,0)
      if(size(co)/=size(u,1)) then
         res%status=LSI_INVALID_INPUT; res%message='ldp: incompatible dimensions'; return
      end if
      if(present(rcond)) then
         call ldp_solve(u,co,res%x,feasible,rcond)
      else
         call ldp_solve(u,co,res%x,feasible)
      end if
      if(.not.feasible) then
         res%status=LSI_INFEASIBLE; res%message='ldp: infeasible constraints'; return
      end if
      res%rnorm=vecnorm(res%x); res%objective=res%rnorm*res%rnorm
   end subroutine ldp

   subroutine ls_ln_vector(a,b,res,rcond,mnorm,x0)
      real(dp), intent(in) :: a(:,:),b(:)
      type(lsi_result), intent(out) :: res
      real(dp), intent(in), optional :: rcond,mnorm(:,:),x0(:)
      integer :: n,st
      n=size(a,2); call init_result(res,n,size(a,1))
      if(size(b)/=size(a,1)) then
         res%status=LSI_INVALID_INPUT; res%message='ls_ln: incompatible dimensions'; return
      end if
      if(present(mnorm)) then
         if(present(x0)) then
            if(present(rcond)) then
               call least_norm_ls(a,b,res%x,res%rank,rcond,mnorm,x0,st)
            else
               call least_norm_ls(a,b,res%x,res%rank,mnorm=mnorm,x0=x0,info=st)
            end if
         else
            if(present(rcond)) then
               call least_norm_ls(a,b,res%x,res%rank,rcond,mnorm=mnorm,info=st)
            else
               call least_norm_ls(a,b,res%x,res%rank,mnorm=mnorm,info=st)
            end if
         end if
      else if(present(x0)) then
         if(present(rcond)) then
            call least_norm_ls(a,b,res%x,res%rank,rcond,x0=x0,info=st)
         else
            call least_norm_ls(a,b,res%x,res%rank,x0=x0,info=st)
         end if
      else
         if(present(rcond)) then
            call least_norm_ls(a,b,res%x,res%rank,rcond,info=st)
         else
            call least_norm_ls(a,b,res%x,res%rank,info=st)
         end if
      end if
      if(st/=0) then
         res%status=LSI_NUMERICAL; res%message='ls_ln: numerical decomposition failed'
      end if
      call finish_ls(a,b,res)
   end subroutine ls_ln_vector

   subroutine ls_ln_multi(a,b,x,rank,status,rcond,mnorm,x0)
      real(dp), intent(in) :: a(:,:),b(:,:)
      real(dp), intent(out) :: x(:,:)
      integer, intent(out) :: rank,status
      real(dp), intent(in), optional :: rcond,mnorm(:,:),x0(:)
      integer :: j,st,rj
      status=0; rank=0; x=0.0_dp
      if(size(b,1)/=size(a,1) .or. size(x,1)/=size(a,2) .or. size(x,2)/=size(b,2)) then
         status=LSI_INVALID_INPUT; return
      end if
      do j=1,size(b,2)
         if(present(mnorm)) then
            if(present(x0)) then
               if(present(rcond)) then
                  call least_norm_ls(a,b(:,j),x(:,j),rj,rcond,mnorm,x0,st)
               else
                  call least_norm_ls(a,b(:,j),x(:,j),rj,mnorm=mnorm,x0=x0,info=st)
               end if
            else
               if(present(rcond)) then
                  call least_norm_ls(a,b(:,j),x(:,j),rj,rcond,mnorm=mnorm,info=st)
               else
                  call least_norm_ls(a,b(:,j),x(:,j),rj,mnorm=mnorm,info=st)
               end if
            end if
         else if(present(x0)) then
            if(present(rcond)) then
               call least_norm_ls(a,b(:,j),x(:,j),rj,rcond,x0=x0,info=st)
            else
               call least_norm_ls(a,b(:,j),x(:,j),rj,x0=x0,info=st)
            end if
         else
            if(present(rcond)) then
               call least_norm_ls(a,b(:,j),x(:,j),rj,rcond,info=st)
            else
               call least_norm_ls(a,b(:,j),x(:,j),rj,info=st)
            end if
         end if
         if(j==1) rank=rj
         if(st/=0) status=LSI_NUMERICAL
      end do
   end subroutine ls_ln_multi

   subroutine ls_ln_svd(a,b,res,rcond)
      real(dp), intent(in) :: a(:,:),b(:)
      type(lsi_result), intent(out) :: res
      real(dp), intent(in), optional :: rcond
      if(present(rcond)) then
         call ls_ln(a,b,res,rcond)
      else
         call ls_ln(a,b,res)
      end if
   end subroutine ls_ln_svd

   subroutine lsi(a,b,res,u,co,rcond)
      real(dp), intent(in) :: a(:,:),b(:)
      type(lsi_result), intent(out) :: res
      real(dp), intent(in), optional :: u(:,:),co(:),rcond
      real(dp), allocatable :: u0(:,:),c0(:),e0(:,:),ce0(:)
      real(dp) :: rc
      integer :: n
      n=size(a,2); rc=1.0e10_dp; if(present(rcond)) rc=rcond
      if(matrix_rank(a,rc)<n) then
         call init_result(res,n,size(a,1)); res%rank=matrix_rank(a,rc)
         res%status=LSI_RANK_DEFICIENT; res%message='lsi: rank deficient design matrix'; return
      end if
      call make_optional_constraints(n,u,co,u0,c0)
      allocate(e0(0,n),ce0(0))
      call constrained_ls(a,b,u0,c0,e0,ce0,res,rc)
   end subroutine lsi

   subroutine lsi_ln(a,b,res,u,co,rcond,mnorm,x0)
      real(dp), intent(in) :: a(:,:),b(:)
      type(lsi_result), intent(out) :: res
      real(dp), intent(in), optional :: u(:,:),co(:),rcond,mnorm(:,:),x0(:)
      real(dp), allocatable :: u0(:,:),c0(:),e0(:,:),ce0(:)
      real(dp) :: rc
      integer :: n
      n=size(a,2); rc=1.0e10_dp; if(present(rcond)) rc=rcond
      call make_optional_constraints(n,u,co,u0,c0)
      if(size(u0,1)==0) then
         if(present(mnorm)) then
            if(present(x0)) then; call ls_ln(a,b,res,rc,mnorm,x0)
            else; call ls_ln(a,b,res,rc,mnorm=mnorm); end if
         else if(present(x0)) then; call ls_ln(a,b,res,rc,x0=x0)
         else; call ls_ln(a,b,res,rc); end if
         return
      end if
      allocate(e0(0,n),ce0(0))
      if(present(mnorm)) then
         if(present(x0)) then
            call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,mnorm,x0)
         else
            call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,mnorm=mnorm)
         end if
      else if(present(x0)) then
         call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,x0=x0)
      else
         call constrained_ls(a,b,u0,c0,e0,ce0,res,rc)
      end if
   end subroutine lsi_ln

   subroutine lsie_ln(a,b,res,u,co,e,ce,rcond,mnorm,x0)
      real(dp), intent(in) :: a(:,:),b(:)
      type(lsi_result), intent(out) :: res
      real(dp), intent(in), optional :: u(:,:),co(:),e(:,:),ce(:),rcond,mnorm(:,:),x0(:)
      real(dp), allocatable :: u0(:,:),c0(:),e0(:,:),ce0(:)
      real(dp) :: rc
      integer :: n
      n=size(a,2); rc=1.0e10_dp; if(present(rcond)) rc=rcond
      call make_optional_constraints(n,u,co,u0,c0)
      call make_optional_equalities(n,e,ce,e0,ce0)
      if(present(mnorm)) then
         if(present(x0)) then
            call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,mnorm,x0)
         else
            call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,mnorm=mnorm)
         end if
      else if(present(x0)) then
         call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,x0=x0)
      else
         call constrained_ls(a,b,u0,c0,e0,ce0,res,rc)
      end if
   end subroutine lsie_ln

   subroutine constrained_ls(a,b,u,co,e,ce,res,rcond,mnorm,x0,lambda_reg)
      real(dp), intent(in) :: a(:,:),b(:),u(:,:),co(:),e(:,:),ce(:),rcond
      type(lsi_result), intent(out) :: res
      real(dp), intent(in), optional :: mnorm(:,:),x0(:),lambda_reg
      real(dp), allocatable :: h(:,:),d(:),mtm(:,:),xbase(:),x(:),z(:,:),ep(:,:),ue(:,:),rhsu(:)
      real(dp), allocatable :: kkt(:,:),rhs(:),sol(:),lami(:)
      logical, allocatable :: active(:)
      real(dp) :: tol,eta,scale,alpha,slack0,slack1,reg
      integer :: n,m,ne,nu,ranke,ranka,st,iter,na,i,j,k,block,drop,maxit
      logical :: feasible
      n=size(a,2); m=size(a,1); ne=size(e,1); nu=size(u,1)
      call init_result(res,n,m); tol=1.0_dp/rcond
      if(size(b)/=m .or. size(u,2)/=n .or. size(co)/=nu .or. size(e,2)/=n .or. size(ce)/=ne) then
         res%status=LSI_INVALID_INPUT; res%message='constrained least squares: incompatible dimensions'; return
      end if
      ranka=matrix_rank(a,rcond); res%rank=ranka
      allocate(h(n,n),d(n),mtm(n,n),xbase(n)); h=matmul(transpose(a),a); d=matmul(transpose(a),b)
      mtm=identity_matrix(n); xbase=0.0_dp
      if(present(mnorm)) mtm=matmul(transpose(mnorm),mnorm)
      if(present(x0)) xbase=x0
      scale=max(1.0_dp,maxval(abs(h)))
      reg=0.0_dp; if(present(lambda_reg)) reg=max(0.0_dp,lambda_reg)
      if(reg>0.0_dp) then
         h=h+(reg*reg)*mtm; d=d+(reg*reg)*matmul(mtm,xbase)
      else if(ranka<n) then
         eta=max(1.0e-13_dp,1000.0_dp*epsilon(1.0_dp))*scale
         h=h+eta*mtm; d=d+eta*matmul(mtm,xbase)
      end if
      allocate(x(n)); x=0.0_dp
      if(ne>0) then
         call equality_feasible(e,ce,x,z,ranke,rcond,st)
         if(st/=0) then
            res%status=LSI_INFEASIBLE; res%message='inconsistent equality constraints'; return
         end if
         if(nu>0 .and. any(matmul(u,x)-co < -100.0_dp*tol)) then
            if(size(z,2)==0) then
               res%status=LSI_INFEASIBLE; res%message='inequalities conflict with equalities'; return
            end if
            allocate(ue(nu,size(z,2)),rhsu(nu),sol(size(z,2)))
            ue=matmul(u,z); rhsu=co-matmul(u,x)
            call ldp_solve(ue,rhsu,sol,feasible,rcond)
            if(.not.feasible) then
               res%status=LSI_INFEASIBLE; res%message='infeasible constraints'; return
            end if
            x=x+matmul(z,sol); deallocate(ue,rhsu,sol)
         end if
      else
         if(nu>0 .and. any(-co < -100.0_dp*tol)) then
            call ldp_solve(u,co,x,feasible,rcond)
            if(.not.feasible) then
               res%status=LSI_INFEASIBLE; res%message='infeasible inequalities'; return
            end if
         end if
      end if
      allocate(active(nu)); active=.false.
      do i=1,nu
         if(abs(dot_product(u(i,:),x)-co(i))<=1000.0_dp*tol*max(1.0_dp,abs(co(i)))) active(i)=.true.
      end do
      maxit=max(100,30*(n+nu+ne+1))
      do iter=1,maxit
         na=count(active); allocate(kkt(n+ne+na,n+ne+na),rhs(n+ne+na),sol(n+ne+na),lami(na))
         kkt=0.0_dp; rhs=0.0_dp; kkt(1:n,1:n)=h; rhs(1:n)=d
         if(ne>0) then
            kkt(1:n,n+1:n+ne)=transpose(e); kkt(n+1:n+ne,1:n)=e; rhs(n+1:n+ne)=ce
         end if
         k=0
         do i=1,nu
            if(active(i)) then
               k=k+1
               kkt(1:n,n+ne+k)=-u(i,:); kkt(n+ne+k,1:n)=u(i,:); rhs(n+ne+k)=co(i)
            end if
         end do
         call dense_solve(kkt,rhs,sol,st)
         if(st/=0) then
            do i=1,n
               kkt(i,i)=kkt(i,i)+1.0e-11_dp*scale
            end do
            call dense_solve(kkt,rhs,sol,st)
         end if
         if(st/=0) then
            res%status=LSI_NUMERICAL; res%message='KKT solve failed'; return
         end if
         if(na>0) lami=sol(n+ne+1:n+ne+na)
         block=0; alpha=1.0_dp
         do i=1,nu
            if(active(i)) cycle
            slack0=dot_product(u(i,:),x)-co(i)
            slack1=dot_product(u(i,:),sol(1:n))-co(i)
            if(slack1 < -100.0_dp*tol .and. slack1<slack0) then
               if(slack0/(slack0-slack1)<alpha) then
                  alpha=max(0.0_dp,slack0/(slack0-slack1)); block=i
               end if
            end if
         end do
         if(block>0 .and. alpha<1.0_dp-100.0_dp*epsilon(1.0_dp)) then
            x=x+alpha*(sol(1:n)-x); active(block)=.true.
            deallocate(kkt,rhs,sol,lami); cycle
         end if
         x=sol(1:n); drop=0; k=0
         do i=1,nu
            if(active(i)) then
               k=k+1
               if(lami(k)<-1000.0_dp*tol) then
                  if(drop==0) then
                     drop=i
                  else
                     j=count(active(1:i-1))
                     if(lami(k)<lami(j)) drop=i
                  end if
               end if
            end if
         end do
         deallocate(kkt,rhs,sol,lami)
         if(drop==0) exit
         active(drop)=.false.
      end do
      if(iter>maxit) then
         res%status=LSI_NUMERICAL; res%message='active-set iteration limit'; return
      end if
      if(ranka<n .and. reg<=0.0_dp) then
         if(present(mnorm)) then
            if(present(x0)) then
               call refine_min_norm_active(a,b,u,co,e,ce,active,x,rcond,mnorm,x0,st)
            else
               call refine_min_norm_active(a,b,u,co,e,ce,active,x,rcond,mnorm=mnorm,status=st)
            end if
         else if(present(x0)) then
            call refine_min_norm_active(a,b,u,co,e,ce,active,x,rcond,x0=x0,status=st)
         else
            call refine_min_norm_active(a,b,u,co,e,ce,active,x,rcond,status=st)
         end if
         if(st/=0) then
            res%status=LSI_NUMERICAL; res%message='least-norm refinement failed'; return
         end if
      end if
      if(ne>0) then
         allocate(ep(n,ne))
         call pseudoinverse(e,ep,ranke,rcond,st)
         x=x+matmul(ep,ce-matmul(e,x))
         deallocate(ep)
      end if
      res%x=x
      where(abs(res%x)<100.0_dp*epsilon(1.0_dp)) res%x=0.0_dp
      if(ne>0) then
         if(maxval(abs(matmul(e,res%x)-ce))>1000.0_dp*tol*max(1.0_dp,maxval(abs(ce)))) then
            res%status=LSI_INFEASIBLE; res%message='equality residual too large'; return
         end if
      end if
      if(nu>0) then
         if(minval(matmul(u,res%x)-co)<-1000.0_dp*tol*max(1.0_dp,maxval(abs(co)))) then
            res%status=LSI_INFEASIBLE; res%message='inequality residual too large'; return
         end if
      end if
      call finish_ls(a,b,res)
   end subroutine constrained_ls

   subroutine refine_min_norm_active(a,b,u,co,e,ce,active,x,rcond,mnorm,x0,status)
      real(dp), intent(in) :: a(:,:),b(:),u(:,:),co(:),e(:,:),ce(:),rcond
      logical, intent(inout) :: active(:)
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in), optional :: mnorm(:,:),x0(:)
      integer, intent(out) :: status
      real(dp), allocatable :: cmat(:,:),crhs(:),xc(:),z(:,:),az(:,:),bz(:),paz(:,:),zls(:)
      real(dp), allocatable :: nz(:,:),dirs(:,:),wmat(:,:),wrhs(:),pwm(:,:),wcoef(:),base(:)
      real(dp), allocatable :: slack(:),pc(:,:),corr(:)
      integer :: n,ne,nu,na,i,k,rankc,rankaz,rankn,rankw,st,bad,loop
      real(dp) :: tol
      n=size(a,2); ne=size(e,1); nu=size(u,1); tol=1.0_dp/rcond; status=0
      do loop=1,nu+1
         na=count(active); allocate(cmat(ne+na,n),crhs(ne+na)); cmat=0.0_dp; crhs=0.0_dp
         if(ne>0) then; cmat(1:ne,:)=e; crhs(1:ne)=ce; end if
         k=ne
         do i=1,nu
            if(active(i)) then; k=k+1; cmat(k,:)=u(i,:); crhs(k)=co(i); end if
         end do
         allocate(xc(n))
         if(size(cmat,1)>0) then
            call equality_feasible(cmat,crhs,xc,z,rankc,rcond,st)
            if(st/=0) then; status=1; return; end if
         else
            xc=0.0_dp; allocate(z(n,n)); z=identity_matrix(n); rankc=0
         end if
         if(size(z,2)==0) then
            x=xc
         else
            allocate(az(size(a,1),size(z,2)),bz(size(b)),paz(size(z,2),size(a,1)),zls(size(z,2)))
            az=matmul(a,z); bz=b-matmul(a,xc)
            call pseudoinverse(az,paz,rankaz,rcond,st); zls=matmul(paz,bz)
            x=xc+matmul(z,zls)
            call null_space(az,nz,rankn,rcond,st)
            if(size(nz,2)>0) then
               dirs=matmul(z,nz); allocate(base(n)); base=x
               if(present(x0)) base=x-x0
               if(present(mnorm)) then
                  allocate(wmat(size(mnorm,1),size(dirs,2)),wrhs(size(mnorm,1)))
                  wmat=matmul(mnorm,dirs); wrhs=-matmul(mnorm,base)
               else
                  allocate(wmat(n,size(dirs,2)),wrhs(n)); wmat=dirs; wrhs=-base
               end if
               allocate(pwm(size(dirs,2),size(wmat,1)),wcoef(size(dirs,2)))
               call pseudoinverse(wmat,pwm,rankw,rcond,st); wcoef=matmul(pwm,wrhs)
               x=x+matmul(dirs,wcoef)
            end if
         end if
         if(size(cmat,1)>0) then
            allocate(pc(n,size(cmat,1)),corr(n))
            call pseudoinverse(cmat,pc,rankc,rcond,st)
            corr=matmul(pc,crhs-matmul(cmat,x)); x=x+corr
            deallocate(pc,corr)
         end if
         if(nu==0) return
         allocate(slack(nu)); slack=matmul(u,x)-co; bad=0
         do i=1,nu
            if(.not.active(i) .and. slack(i)<-1000.0_dp*tol*max(1.0_dp,abs(co(i)))) then
               if(bad==0 .or. slack(i)<slack(bad)) bad=i
            end if
         end do
         if(bad==0) return
         active(bad)=.true.
         deallocate(cmat,crhs,xc,z,slack)
         if(allocated(az)) deallocate(az,bz,paz,zls)
         if(allocated(nz)) deallocate(nz)
         if(allocated(dirs)) deallocate(dirs)
         if(allocated(base)) deallocate(base)
         if(allocated(wmat)) deallocate(wmat,wrhs,pwm,wcoef)
      end do
      status=1
   end subroutine refine_min_norm_active

   subroutine equality_feasible(e,ce,x,z,ranke,rcond,status)
      real(dp), intent(in) :: e(:,:),ce(:),rcond
      real(dp), intent(out) :: x(:)
      real(dp), allocatable, intent(out) :: z(:,:)
      integer, intent(out) :: ranke,status
      real(dp), allocatable :: ep(:,:)
      integer :: st
      status=0; allocate(ep(size(e,2),size(e,1)))
      call pseudoinverse(e,ep,ranke,rcond,st); x=matmul(ep,ce)
      if(maxval(abs(matmul(e,x)-ce))>100.0_dp/rcond*max(1.0_dp,maxval(abs(ce)))) then
         status=1; allocate(z(size(e,2),0)); return
      end if
      call null_space(e,z,ranke,rcond,st)
   end subroutine equality_feasible

   subroutine lsi_reg(a,b,res,u,co,rcond,mnorm,x0)
      real(dp), intent(in) :: a(:,:),b(:)
      type(lsi_result), intent(out) :: res
      real(dp), intent(in), optional :: u(:,:),co(:),rcond,mnorm(:,:),x0(:)
      real(dp), allocatable :: ata(:,:),ev(:),vec(:,:),u0(:,:),c0(:),e0(:,:),ce0(:)
      real(dp) :: rc,smax,lambda
      integer :: n,rank,st
      n=size(a,2); rc=1.0e10_dp; if(present(rcond)) rc=rcond
      rank=matrix_rank(a,rc); allocate(ata(n,n),ev(n),vec(n,n)); ata=matmul(transpose(a),a)
      call symmetric_eigen(ata,ev,vec,st); smax=sqrt(max(0.0_dp,ev(1))); lambda=0.0_dp
      call make_optional_constraints(n,u,co,u0,c0); allocate(e0(0,n),ce0(0))
      if(rank>0 .and. rank<n) lambda=smax/sqrt(rc)
      if(lambda<=0.0_dp) then
         if(size(u0,1)>0) then
            if(present(mnorm)) then
               if(present(x0)) then; call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,mnorm,x0)
               else; call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,mnorm=mnorm); end if
            else if(present(x0)) then; call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,x0=x0)
            else; call constrained_ls(a,b,u0,c0,e0,ce0,res,rc); end if
         else
            if(present(mnorm)) then
               if(present(x0)) then; call ls_ln(a,b,res,rc,mnorm,x0)
               else; call ls_ln(a,b,res,rc,mnorm=mnorm); end if
            else if(present(x0)) then; call ls_ln(a,b,res,rc,x0=x0)
            else; call ls_ln(a,b,res,rc); end if
         end if
         res%lambda=0.0_dp; return
      end if
      if(present(mnorm)) then
         if(present(x0)) then
            call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,mnorm,x0,lambda)
         else
            call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,mnorm=mnorm,lambda_reg=lambda)
         end if
      else if(present(x0)) then
         call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,x0=x0,lambda_reg=lambda)
      else
         call constrained_ls(a,b,u0,c0,e0,ce0,res,rc,lambda_reg=lambda)
      end if
      res%lambda=lambda
   end subroutine lsi_reg

   subroutine tls(a,b,x,status)
      real(dp), intent(in) :: a(:,:),b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: ab(:,:),ata(:,:),ev(:),v(:,:)
      integer :: n,st
      n=size(a,2); x=0.0_dp; status=0
      if(size(b)/=size(a,1) .or. size(x)/=n) then; status=1; return; end if
      allocate(ab(size(a,1),n+1),ata(n+1,n+1),ev(n+1),v(n+1,n+1))
      ab(:,1:n)=a; ab(:,n+1)=b; ata=matmul(transpose(ab),ab)
      call symmetric_eigen(ata,ev,v,st)
      if(abs(v(n+1,n+1))<=100.0_dp*epsilon(1.0_dp)) then; status=2; return; end if
      x=v(1:n,n+1)/(-v(n+1,n+1))
   end subroutine tls

   subroutine nulla_matrix(m,basis,rank,rcond)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable, intent(out) :: basis(:,:)
      integer, intent(out) :: rank
      real(dp), intent(in), optional :: rcond
      integer :: st
      if(present(rcond)) then
         call null_space(transpose(m),basis,rank,rcond,st)
      else
         call null_space(transpose(m),basis,rank,info=st)
      end if
   end subroutine nulla_matrix

   subroutine nulla_vector(v,basis,rank,rcond)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable, intent(out) :: basis(:,:)
      integer, intent(out) :: rank
      real(dp), intent(in), optional :: rcond
      real(dp), allocatable :: m(:,:)
      allocate(m(size(v),1)); m(:,1)=v
      if(present(rcond)) then
         call nulla_matrix(m,basis,rank,rcond)
      else
         call nulla_matrix(m,basis,rank)
      end if
   end subroutine nulla_vector

   subroutine pnull(a,b,xp,basis,rank,rcond,status)
      real(dp), intent(in) :: a(:,:),b(:)
      real(dp), intent(out) :: xp(:)
      real(dp), allocatable, intent(out) :: basis(:,:)
      integer, intent(out) :: rank,status
      real(dp), intent(in), optional :: rcond
      real(dp), allocatable :: ap(:,:)
      integer :: st
      status=0; allocate(ap(size(a,2),size(a,1)))
      if(present(rcond)) then
         call pseudoinverse(a,ap,rank,rcond,st); call null_space(a,basis,rank,rcond,st)
      else
         call pseudoinverse(a,ap,rank,info=st); call null_space(a,basis,rank,info=st)
      end if
      xp=matmul(ap,b); if(st/=0) status=st
   end subroutine pnull

   subroutine uplo_to_uco(lower,upper,u,co)
      real(dp), intent(in), optional :: lower(:),upper(:)
      real(dp), allocatable, intent(out) :: u(:,:),co(:)
      integer :: n,nl,nu,i,k
      n=0; if(present(lower)) n=size(lower); if(present(upper)) n=max(n,size(upper))
      nl=0; nu=0
      if(present(lower)) nl=count(ieee_is_finite(lower))
      if(present(upper)) nu=count(ieee_is_finite(upper))
      allocate(u(nl+nu,n),co(nl+nu)); u=0.0_dp; co=0.0_dp; k=0
      if(present(upper)) then
         do i=1,n
            if(ieee_is_finite(upper(i))) then; k=k+1; u(k,i)=-1.0_dp; co(k)=-upper(i); end if
         end do
      end if
      if(present(lower)) then
         do i=1,n
            if(ieee_is_finite(lower(i))) then; k=k+1; u(k,i)=1.0_dp; co(k)=lower(i); end if
         end do
      end if
   end subroutine uplo_to_uco

   subroutine make_optional_constraints(n,u,co,u0,c0)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: u(:,:),co(:)
      real(dp), allocatable, intent(out) :: u0(:,:),c0(:)
      if(present(u)) then
         allocate(u0(size(u,1),n)); u0=u
         if(present(co)) then; allocate(c0(size(co))); c0=co
         else; allocate(c0(size(u,1))); c0=0.0_dp; end if
      else
         allocate(u0(0,n),c0(0))
      end if
   end subroutine make_optional_constraints

   subroutine make_optional_equalities(n,e,ce,e0,ce0)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: e(:,:),ce(:)
      real(dp), allocatable, intent(out) :: e0(:,:),ce0(:)
      if(present(e)) then
         allocate(e0(size(e,1),n)); e0=e
         if(present(ce)) then; allocate(ce0(size(ce))); ce0=ce
         else; allocate(ce0(size(e,1))); ce0=0.0_dp; end if
      else
         allocate(e0(0,n),ce0(0))
      end if
   end subroutine make_optional_equalities
end module nlsic_linear
