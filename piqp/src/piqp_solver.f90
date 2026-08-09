! SPDX-License-Identifier: BSD-2-Clause
module piqp_solver
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf
   use piqp_kinds, only : dp
   use piqp_types, only : piqp_settings_type, piqp_result_type, piqp_info_type, PIQP_SOLVED, &
      PIQP_MAX_ITER_REACHED, PIQP_PRIMAL_INFEASIBLE, PIQP_NUMERICS, &
      PIQP_INVALID_SETTINGS, PIQP_UNSOLVED
   use piqp_linalg, only : chol_solve_spd, norm_inf, all_finite
   implicit none
   private

   type, public :: piqp_model_type
      real(dp), allocatable :: pmat(:,:), c(:), amat(:,:), b(:), gmat(:,:)
      real(dp), allocatable :: h_l(:), h_u(:), x_l(:), x_u(:)
      type(piqp_settings_type) :: settings
      type(piqp_result_type) :: result
      logical :: setup_done = .false.
   contains
      procedure :: setup => piqp_model_setup
      procedure :: solve => piqp_model_solve
      procedure :: update => piqp_model_update
      procedure :: get_dims => piqp_model_get_dims
   end type piqp_model_type

   public :: solve_piqp_dense

contains

   subroutine solve_piqp_dense(pmat, c, result, amat, b, gmat, h_l, h_u, x_l, x_u, settings)
      real(dp), intent(in), optional :: pmat(:,:)
      real(dp), intent(in) :: c(:)
      type(piqp_result_type), intent(out) :: result
      real(dp), intent(in), optional :: amat(:,:), b(:), gmat(:,:), h_l(:), h_u(:), x_l(:), x_u(:)
      type(piqp_settings_type), intent(in), optional :: settings
      type(piqp_model_type) :: model
      integer :: n
      n = size(c)
      if (present(pmat)) then
         call model%setup(pmat=pmat, c=c, amat=amat, b=b, gmat=gmat, h_l=h_l, h_u=h_u, &
            x_l=x_l, x_u=x_u, settings=settings)
      else
         call model%setup(pmat=reshape([real(dp)::], [0,0]), c=c, amat=amat, b=b, gmat=gmat, &
            h_l=h_l, h_u=h_u, x_l=x_l, x_u=x_u, settings=settings, zero_p_n=n)
      end if
      call model%solve()
      result = model%result
   end subroutine solve_piqp_dense

   subroutine piqp_model_setup(self, pmat, c, amat, b, gmat, h_l, h_u, x_l, x_u, settings, zero_p_n)
      class(piqp_model_type), intent(inout) :: self
      real(dp), intent(in) :: pmat(:,:), c(:)
      real(dp), intent(in), optional :: amat(:,:), b(:), gmat(:,:), h_l(:), h_u(:), x_l(:), x_u(:)
      type(piqp_settings_type), intent(in), optional :: settings
      integer, intent(in), optional :: zero_p_n
      integer :: n, p, m
      real(dp) :: inf
      inf = ieee_value(0.0_dp, ieee_positive_inf)
      n = size(c)
      if (present(zero_p_n)) then
         if (zero_p_n /= n) then
            self%setup_done = .false.
            return
         end if
         allocate(self%pmat(n,n), source=0.0_dp)
      else
         if (size(pmat,1) /= n .or. size(pmat,2) /= n) then
            self%setup_done = .false.
            return
         end if
         allocate(self%pmat(n,n))
         self%pmat = 0.5_dp*(pmat + transpose(pmat))
      end if
      allocate(self%c(n)); self%c = c
      p = 0
      if (present(amat)) p = size(amat,1)
      allocate(self%amat(p,n), self%b(p))
      self%amat = 0.0_dp; self%b = 0.0_dp
      if (p > 0) then
         if (size(amat,2) /= n .or. .not. present(b) .or. size(b) /= p) then
            self%setup_done = .false.; return
         end if
         self%amat = amat; self%b = b
      end if
      m = 0
      if (present(gmat)) m = size(gmat,1)
      allocate(self%gmat(m,n), self%h_l(m), self%h_u(m))
      self%gmat = 0.0_dp; self%h_l = -inf; self%h_u = inf
      if (m > 0) then
         if (size(gmat,2) /= n) then
            self%setup_done = .false.; return
         end if
         self%gmat = gmat
         if (present(h_l)) then
            if (size(h_l) /= m) then; self%setup_done=.false.; return; end if
            self%h_l = h_l
         end if
         if (present(h_u)) then
            if (size(h_u) /= m) then; self%setup_done=.false.; return; end if
            self%h_u = h_u
         end if
      end if
      allocate(self%x_l(n), self%x_u(n))
      self%x_l = -inf; self%x_u = inf
      if (present(x_l)) then
         if (size(x_l) /= n) then; self%setup_done=.false.; return; end if
         self%x_l = x_l
      end if
      if (present(x_u)) then
         if (size(x_u) /= n) then; self%setup_done=.false.; return; end if
         self%x_u = x_u
      end if
      if (present(settings)) self%settings = settings
      self%setup_done = .true.
      self%result%info%status = -9
   end subroutine piqp_model_setup

   subroutine piqp_model_update(self, pmat, c, amat, b, gmat, h_l, h_u, x_l, x_u, settings, info)
      class(piqp_model_type), intent(inout) :: self
      real(dp), intent(in), optional :: pmat(:,:), c(:), amat(:,:), b(:), gmat(:,:), h_l(:), h_u(:), x_l(:), x_u(:)
      type(piqp_settings_type), intent(in), optional :: settings
      integer, intent(out), optional :: info
      integer :: n, p, m, ierr
      ierr = 0
      if (.not. self%setup_done) then
         ierr = 1
         if (present(info)) info = ierr
         return
      end if
      n = size(self%c); p = size(self%b); m = size(self%h_l)
      if (present(pmat)) then
         if (size(pmat,1)/=n .or. size(pmat,2)/=n) ierr=2
         if (ierr==0) self%pmat = 0.5_dp*(pmat+transpose(pmat))
      end if
      if (present(c) .and. ierr==0) then
         if (size(c)/=n) ierr=2
         if (ierr==0) self%c=c
      end if
      if (present(amat) .and. ierr==0) then
         if (size(amat,1)/=p .or. size(amat,2)/=n) ierr=2
         if (ierr==0) self%amat=amat
      end if
      if (present(b) .and. ierr==0) then
         if (size(b)/=p) ierr=2
         if (ierr==0) self%b=b
      end if
      if (present(gmat) .and. ierr==0) then
         if (size(gmat,1)/=m .or. size(gmat,2)/=n) ierr=2
         if (ierr==0) self%gmat=gmat
      end if
      if (present(h_l) .and. ierr==0) then
         if (size(h_l)/=m) ierr=2
         if (ierr==0) self%h_l=h_l
      end if
      if (present(h_u) .and. ierr==0) then
         if (size(h_u)/=m) ierr=2
         if (ierr==0) self%h_u=h_u
      end if
      if (present(x_l) .and. ierr==0) then
         if (size(x_l)/=n) ierr=2
         if (ierr==0) self%x_l=x_l
      end if
      if (present(x_u) .and. ierr==0) then
         if (size(x_u)/=n) ierr=2
         if (ierr==0) self%x_u=x_u
      end if
      if (present(settings) .and. ierr==0) self%settings=settings
      if (present(info)) info=ierr
   end subroutine piqp_model_update

   subroutine piqp_model_get_dims(self, n, p, m)
      class(piqp_model_type), intent(in) :: self
      integer, intent(out) :: n, p, m
      if (.not. self%setup_done) then
         n=0; p=0; m=0
      else
         n=size(self%c); p=size(self%b); m=size(self%h_l)
      end if
   end subroutine piqp_model_get_dims

   subroutine piqp_model_solve(self)
      class(piqp_model_type), intent(inout) :: self
      if (.not. self%setup_done) then
         self%result%info%status = -9
         return
      end if
      call solve_core(self%pmat, self%c, self%amat, self%b, self%gmat, self%h_l, self%h_u, &
         self%x_l, self%x_u, self%settings, self%result)
   end subroutine piqp_model_solve

   subroutine solve_core(pmat, cvec, amat, bvec, gmat, h_l, h_u, x_l, x_u, settings, result)
      real(dp), intent(in) :: pmat(:,:), cvec(:), amat(:,:), bvec(:), gmat(:,:), h_l(:), h_u(:), x_l(:), x_u(:)
      type(piqp_settings_type), intent(in) :: settings
      type(piqp_result_type), intent(inout) :: result
      real(dp), allocatable :: cmat(:,:), d(:), x(:), y(:), s(:), z(:), xi(:), lam(:)
      real(dp), allocatable :: rdual(:), req(:), rineq(:), dx(:), dy(:), ds(:), dz(:)
      real(dp), allocatable :: dsa(:), dza(:), rhscent(:), work(:)
      integer, allocatable :: ckind(:), cidx(:)
      integer :: n, p, m, q, iter, linfo, i
      real(dp) :: rho, delta, mu, mu_aff, sigma, ap, ad, pstep, dstep, tau
      real(dp) :: pres, dres, gap, pobj, dobj, pden, dden, oldpres, olddres, muold, murate
      real(dp) :: reg_limit, t0, t1
      logical :: solved
      n=size(cvec); p=size(bvec); m=size(h_l)
      call allocate_result(result,n,p,m)
      result%info%status = PIQP_UNSOLVED
      if (.not. settings%valid()) then
         result%info%status = PIQP_INVALID_SETTINGS
         return
      end if
      if (any((ieee_is_finite(x_l) .and. ieee_is_finite(x_u)) .and. x_l > x_u) .or. &
          any((ieee_is_finite(h_l) .and. ieee_is_finite(h_u)) .and. h_l > h_u)) then
         result%info%status = PIQP_PRIMAL_INFEASIBLE
         return
      end if
      call build_inequalities(gmat,h_l,h_u,x_l,x_u,cmat,d,ckind,cidx)
      q=size(d)
      allocate(x(n),y(p),s(q),z(q),xi(n),lam(p),rdual(n),req(p),rineq(q),dx(n),dy(p),ds(q),dz(q), &
         dsa(q),dza(q),rhscent(q),work(max(1,q)))
      rho=settings%rho_init; delta=settings%delta_init; reg_limit=settings%reg_lower_limit; tau=settings%tau
      call cpu_time(t0)
      call initial_point(pmat,cvec,amat,bvec,rho,delta,settings,x,y,linfo)
      if (linfo /= 0) then
         result%info%status=PIQP_NUMERICS; return
      end if
      xi=x; lam=y
      if (q>0) then
         s = d - matmul(cmat,x)
         do i=1,q
            s(i)=max(1.0_dp,s(i))
         end do
         z=1.0_dp
      end if
      oldpres=huge(1.0_dp); olddres=huge(1.0_dp)
      result%info%rho=rho; result%info%delta=delta; result%info%reg_limit=reg_limit
      solved=.false.
      do iter=0,settings%max_iter
         call residual_metrics(pmat,cvec,amat,bvec,cmat,d,x,y,z,pres,dres,gap,pobj,dobj,pden,dden)
         result%info%prev_primal_res=result%info%primal_res
         result%info%prev_dual_res=result%info%dual_res
         result%info%primal_res=pres
         result%info%dual_res=dres
         result%info%primal_res_rel=pres/max(1.0_dp,pden)
         result%info%dual_res_rel=dres/max(1.0_dp,dden)
         result%info%primal_obj=pobj; result%info%dual_obj=dobj
         result%info%duality_gap=gap
         result%info%duality_gap_rel=gap/max(1.0_dp,abs(pobj),abs(dobj))
         result%info%iter=iter
         if (q>0) then
            mu=dot_product(s,z)/real(q,dp)
         else
            mu=0.0_dp
         end if
         result%info%mu=mu
         if ((pres < settings%eps_abs .or. result%info%primal_res_rel < settings%eps_rel) .and. &
             (dres < settings%eps_abs .or. result%info%dual_res_rel < settings%eps_rel) .and. &
             (.not.settings%check_duality_gap .or. gap < settings%eps_duality_gap_abs .or. &
              result%info%duality_gap_rel < settings%eps_duality_gap_rel)) then
            solved=.true.; exit
         end if
         if (iter == settings%max_iter) exit
         rdual=matmul(pmat,x)+cvec
         if (p>0) rdual=rdual+matmul(transpose(amat),y)
         if (q>0) rdual=rdual+matmul(transpose(cmat),z)
         rdual=rdual+rho*(x-xi)
         if (p>0) req=matmul(amat,x)-bvec-delta*(y-lam)
         if (q>0) rineq=matmul(cmat,x)+s-d
         if (q>0) then
            rhscent=-s*z
            call newton_step(pmat,amat,cmat,s,z,rdual,req,rineq,rhscent,rho,delta,settings,dx,dy,ds,dz,linfo)
            if (linfo/=0) then
               result%info%status=PIQP_NUMERICS; return
            end if
            call max_step(s,ds,ap); call max_step(z,dz,ad)
            mu_aff=dot_product(s+ap*ds,z+ad*dz)/real(q,dp)
            if (mu>0.0_dp) then
               sigma=max(0.0_dp,min(1.0_dp,mu_aff/mu))**3
            else
               sigma=0.0_dp
            end if
            dsa=ds; dza=dz
            rhscent=-s*z-dsa*dza+sigma*mu
            call newton_step(pmat,amat,cmat,s,z,rdual,req,rineq,rhscent,rho,delta,settings,dx,dy,ds,dz,linfo)
            if (linfo/=0) then
               result%info%status=PIQP_NUMERICS; return
            end if
            call max_step(s,ds,ap); call max_step(z,dz,ad)
            pstep=min(1.0_dp,tau*ap); dstep=min(1.0_dp,tau*ad)
            x=x+pstep*dx; if(p>0) y=y+dstep*dy
            s=s+pstep*ds; z=z+dstep*dz
            result%info%sigma=sigma; result%info%primal_step=pstep; result%info%dual_step=dstep
            muold=max(mu,tiny(1.0_dp)); mu=dot_product(s,z)/real(q,dp)
            murate=max(0.0_dp,(muold-mu)/muold)
         else
            rhscent=[real(dp)::]
            call newton_step(pmat,amat,cmat,s,z,rdual,req,rineq,rhscent,rho,delta,settings,dx,dy,ds,dz,linfo)
            if (linfo/=0) then; result%info%status=PIQP_NUMERICS; return; end if
            x=x+dx; if(p>0) y=y+dy
            result%info%primal_step=1.0_dp; result%info%dual_step=1.0_dp; murate=0.5_dp
         end if
         if (.not.all_finite(x) .or. .not.all_finite(y) .or. .not.all_finite(s) .or. .not.all_finite(z)) then
            result%info%status=PIQP_NUMERICS; return
         end if
         call residual_metrics(pmat,cvec,amat,bvec,cmat,d,x,y,z,pres,dres,gap,pobj,dobj,pden,dden)
         if (dres < 0.95_dp*olddres .or. dres < settings%eps_abs) then
            xi=x; rho=max(reg_limit,(1.0_dp-murate)*rho); result%info%no_primal_update=0
         else
            result%info%no_primal_update=result%info%no_primal_update+1
            rho=max(reg_limit,(1.0_dp-0.666_dp*murate)*rho)
         end if
         if (pres < 0.95_dp*oldpres .or. pres < settings%eps_abs) then
            lam=y; delta=max(reg_limit,(1.0_dp-murate)*delta); result%info%no_dual_update=0
         else
            result%info%no_dual_update=result%info%no_dual_update+1
            delta=max(reg_limit,(1.0_dp-0.666_dp*murate)*delta)
         end if
         oldpres=pres; olddres=dres
         result%info%rho=rho; result%info%delta=delta
         if (settings%verbose) then
            write(*,'(i4,2x,es13.5,2x,es10.3,2x,es10.3,2x,es10.3)') iter,pobj,pres,dres,mu
         end if
      end do
      if (solved) then
         result%info%status=PIQP_SOLVED
      else
         result%info%status=PIQP_MAX_ITER_REACHED
      end if
      call unpack_result(x,y,s,z,ckind,cidx,result)
      call cpu_time(t1); result%info%solve_time=t1-t0; result%info%run_time=result%info%solve_time
   end subroutine solve_core

   subroutine initial_point(pmat,cvec,amat,bvec,rho,delta,settings,x,y,info)
      real(dp),intent(in)::pmat(:,:),cvec(:),amat(:,:),bvec(:),rho,delta
      type(piqp_settings_type),intent(in)::settings
      real(dp),intent(out)::x(:),y(:)
      integer,intent(out)::info
      real(dp),allocatable::h(:,:),rhs(:)
      integer::n,p,i
      n=size(cvec); p=size(bvec); allocate(h(n,n),rhs(n))
      h=pmat; do i=1,n; h(i,i)=h(i,i)+rho; end do
      rhs=-cvec
      if(p>0) then
         h=h+matmul(transpose(amat),amat)/delta
         rhs=rhs+matmul(transpose(amat),bvec)/delta
      end if
      call chol_solve_spd(h,rhs,x,info,settings%iterative_refinement_max_iter, &
         settings%iterative_refinement_eps_abs,settings%iterative_refinement_eps_rel)
      if(info==0 .and. p>0) y=(matmul(amat,x)-bvec)/delta
   end subroutine initial_point

   subroutine newton_step(pmat,amat,cmat,s,z,rdual,req,rineq,rhscent,rho,delta,settings,dx,dy,ds,dz,info)
      real(dp),intent(in)::pmat(:,:),amat(:,:),cmat(:,:),s(:),z(:),rdual(:),req(:),rineq(:),rhscent(:),rho,delta
      type(piqp_settings_type),intent(in)::settings
      real(dp),intent(out)::dx(:),dy(:),ds(:),dz(:)
      integer,intent(out)::info
      real(dp),allocatable::h(:,:),rhs(:),w(:),v(:)
      integer::n,p,q,i
      n=size(rdual); p=size(req); q=size(s)
      allocate(h(n,n),rhs(n),w(q),v(q)); h=pmat
      do i=1,n; h(i,i)=h(i,i)+rho; end do
      rhs=-rdual
      if(q>0) then
         w=z/s
         do i=1,q
            h=h+w(i)*outer(cmat(i,:),cmat(i,:))
            v(i)=(rhscent(i)+z(i)*rineq(i))/s(i)
         end do
         rhs=rhs-matmul(transpose(cmat),v)
      end if
      if(p>0) then
         h=h+matmul(transpose(amat),amat)/delta
         rhs=rhs-matmul(transpose(amat),req)/delta
      end if
      call chol_solve_spd(h,rhs,dx,info,settings%iterative_refinement_max_iter, &
         settings%iterative_refinement_eps_abs,settings%iterative_refinement_eps_rel)
      if(info/=0) return
      if(p>0) dy=(matmul(amat,dx)+req)/delta
      if(q>0) then
         ds=-rineq-matmul(cmat,dx)
         dz=(rhscent-z*ds)/s
      end if
   end subroutine newton_step

   function outer(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::i
      do i=1,size(a); c(i,:)=a(i)*b; end do
   end function outer

   subroutine max_step(v,dv,alpha)
      real(dp),intent(in)::v(:),dv(:)
      real(dp),intent(out)::alpha
      integer::i
      alpha=1.0_dp
      do i=1,size(v)
         if(dv(i)<0.0_dp) alpha=min(alpha,-v(i)/dv(i))
      end do
   end subroutine max_step

   subroutine residual_metrics(pmat,cvec,amat,bvec,cmat,d,x,y,z,pres,dres,gap,pobj,dobj,pden,dden)
      real(dp),intent(in)::pmat(:,:),cvec(:),amat(:,:),bvec(:),cmat(:,:),d(:),x(:),y(:),z(:)
      real(dp),intent(out)::pres,dres,gap,pobj,dobj,pden,dden
      real(dp),allocatable::rd(:),re(:),ri(:),cx(:)
      integer::p,q
      p=size(bvec); q=size(d); allocate(rd(size(x)),re(p),ri(q),cx(q))
      rd=matmul(pmat,x)+cvec
      if(p>0) rd=rd+matmul(transpose(amat),y)
      if(q>0) rd=rd+matmul(transpose(cmat),z)
      if(p>0) re=matmul(amat,x)-bvec
      if(q>0) then
         cx=matmul(cmat,x); ri=max(cx-d,0.0_dp)
      end if
      pres=max(norm_inf(re),norm_inf(ri)); dres=norm_inf(rd)
      pden=max(1.0_dp,norm_inf(bvec),norm_inf(d))
      dden=max(1.0_dp,norm_inf(cvec),norm_inf(matmul(pmat,x)))
      pobj=0.5_dp*dot_product(x,matmul(pmat,x))+dot_product(cvec,x)
      dobj=-0.5_dp*dot_product(x,matmul(pmat,x))
      if(p>0) dobj=dobj-dot_product(bvec,y)
      if(q>0) dobj=dobj-dot_product(d,z)
      gap=abs(pobj-dobj)
   end subroutine residual_metrics

   subroutine build_inequalities(gmat,h_l,h_u,x_l,x_u,cmat,d,kind,idx)
      real(dp),intent(in)::gmat(:,:),h_l(:),h_u(:),x_l(:),x_u(:)
      real(dp),allocatable,intent(out)::cmat(:,:),d(:)
      integer,allocatable,intent(out)::kind(:),idx(:)
      integer::m,n,q,i,k
      m=size(h_l); n=size(x_l); q=0
      do i=1,m
         if(ieee_is_finite(h_l(i))) q=q+1
         if(ieee_is_finite(h_u(i))) q=q+1
      end do
      do i=1,n
         if(ieee_is_finite(x_l(i))) q=q+1
         if(ieee_is_finite(x_u(i))) q=q+1
      end do
      allocate(cmat(q,n),d(q),kind(q),idx(q)); cmat=0.0_dp; k=0
      do i=1,m
         if(ieee_is_finite(h_l(i))) then
            k=k+1; cmat(k,:)=-gmat(i,:); d(k)=-h_l(i); kind(k)=2; idx(k)=i
         end if
         if(ieee_is_finite(h_u(i))) then
            k=k+1; cmat(k,:)=gmat(i,:); d(k)=h_u(i); kind(k)=1; idx(k)=i
         end if
      end do
      do i=1,n
         if(ieee_is_finite(x_l(i))) then
            k=k+1; cmat(k,i)=-1.0_dp; d(k)=-x_l(i); kind(k)=4; idx(k)=i
         end if
         if(ieee_is_finite(x_u(i))) then
            k=k+1; cmat(k,i)=1.0_dp; d(k)=x_u(i); kind(k)=3; idx(k)=i
         end if
      end do
   end subroutine build_inequalities

   subroutine allocate_result(r,n,p,m)
      type(piqp_result_type),intent(inout)::r
      integer,intent(in)::n,p,m
      if (allocated(r%x)) deallocate(r%x)
      if (allocated(r%y)) deallocate(r%y)
      if (allocated(r%z_l)) deallocate(r%z_l)
      if (allocated(r%z_u)) deallocate(r%z_u)
      if (allocated(r%z_bl)) deallocate(r%z_bl)
      if (allocated(r%z_bu)) deallocate(r%z_bu)
      if (allocated(r%s_l)) deallocate(r%s_l)
      if (allocated(r%s_u)) deallocate(r%s_u)
      if (allocated(r%s_bl)) deallocate(r%s_bl)
      if (allocated(r%s_bu)) deallocate(r%s_bu)
      r%info = piqp_info_type()
      allocate(r%x(n),r%y(p),r%z_l(m),r%z_u(m),r%z_bl(n),r%z_bu(n), &
         r%s_l(m),r%s_u(m),r%s_bl(n),r%s_bu(n))
      r%x=0.0_dp;r%y=0.0_dp;r%z_l=0.0_dp;r%z_u=0.0_dp;r%z_bl=0.0_dp;r%z_bu=0.0_dp
      r%s_l=0.0_dp;r%s_u=0.0_dp;r%s_bl=0.0_dp;r%s_bu=0.0_dp
   end subroutine allocate_result

   subroutine unpack_result(x,y,s,z,kind,idx,r)
      real(dp),intent(in)::x(:),y(:),s(:),z(:)
      integer,intent(in)::kind(:),idx(:)
      type(piqp_result_type),intent(inout)::r
      real(dp)::inf
      integer::k
      inf=ieee_value(0.0_dp,ieee_positive_inf)
      r%x=x;r%y=y;r%z_l=0.0_dp;r%z_u=0.0_dp;r%z_bl=0.0_dp;r%z_bu=0.0_dp
      r%s_l=inf;r%s_u=inf;r%s_bl=inf;r%s_bu=inf
      do k=1,size(z)
         select case(kind(k))
         case(1);r%z_u(idx(k))=z(k);r%s_u(idx(k))=s(k)
         case(2);r%z_l(idx(k))=z(k);r%s_l(idx(k))=s(k)
         case(3);r%z_bu(idx(k))=z(k);r%s_bu(idx(k))=s(k)
         case(4);r%z_bl(idx(k))=z(k);r%s_bl(idx(k))=s(k)
         end select
      end do
   end subroutine unpack_result
end module piqp_solver
