! SPDX-License-Identifier: GPL-3.0-only
module sgt_other_learning
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use sgt_kinds, only : dp
   use sgt_status, only : sgt_ok, sgt_invalid_input, sgt_no_convergence
   use sgt_types, only : graph_result
   use sgt_linalg, only : symmetric_pseudoinverse, inverse_matrix, solve_linear_system, &
      symmetric_eigen_jacobi, inverse_symmetric_positive_definite, frobenius_norm, &
      vector_norm2, identity_matrix, project_simplex, project_simplex_sum, &
      nnqp_projected_gradient, helmert_basis, rows_correlation, outer_product, &
      trace_matrix, is_positive_definite
   use sgt_operators, only : L, A, D, Linv, Ainv, Lstar, Mmat
   use sgt_utils, only : pairwise_matrix_rownorm2, upper_view_vec
   use sgt_updates, only : initialize_weights_naive, spectral_reconstruction
   use sgt_objectives, only : vanilla_objective
   use sgt_initial_graph, only : build_initial_graph, laplacian_from_directed_affinity
   implicit none
   private
   public :: learn_smooth_approx_graph, cluster_k_component_graph
   public :: learn_smooth_graph, learn_graph_sigrep
   public :: learn_laplacian_gle_mm, learn_laplacian_gle_admm
   public :: learn_combinatorial_graph_laplacian
contains
   pure function adjacency_from_laplacian(lw) result(aw)
      real(dp), intent(in) :: lw(:,:)
      real(dp), allocatable :: aw(:,:)
      integer :: n,i
      n=size(lw,1); allocate(aw(n,n)); aw=-lw
      do i=1,n
         aw(i,i)=0.0_dp
      end do
   end function adjacency_from_laplacian


   subroutine learn_smooth_approx_graph(y,m,result)
      real(dp), intent(in) :: y(:,:)
      integer, intent(in) :: m
      type(graph_result), intent(out) :: result
      real(dp), allocatable :: affinity(:,:)
      integer :: status
      result=graph_result()
      call build_initial_graph(y,m,affinity,status)
      if (status/=sgt_ok) then; result%status=status; return; end if
      result%laplacian=laplacian_from_directed_affinity(affinity)
      result%adjacency=adjacency_from_laplacian(result%laplacian)
      result%convergence=.true.; result%iterations=1; result%status=sgt_ok
   end subroutine learn_smooth_approx_graph

   subroutine cluster_k_component_graph(y,result,k,m,lambda,eigtol,edgetol,maxiter)
      real(dp), intent(in) :: y(:,:)
      type(graph_result), intent(out) :: result
      integer, intent(in), optional :: k,m,maxiter
      real(dp), intent(in), optional :: lambda,eigtol,edgetol
      real(dp), allocatable :: affinity(:,:),s(:,:),ls(:,:),la(:,:),f(:,:),v(:,:),eig(:),eigvec(:,:)
      real(dp), allocatable :: pvec(:),row(:),lseq(:),times(:)
      real(dp) :: lmd,eig_tol,edge_tol,start,now
      integer :: kval,mval,maxit,n,i,iter,status,nzero,nused
      logical :: converged
      result=graph_result()
      kval=1; if (present(k)) kval=k
      mval=5; if (present(m)) mval=m
      maxit=1000; if (present(maxiter)) maxit=max(1,maxiter)
      lmd=1.0_dp; if (present(lambda)) lmd=lambda
      eig_tol=1e-9_dp; if (present(eigtol)) eig_tol=eigtol
      edge_tol=1e-6_dp; if (present(edgetol)) edge_tol=edgetol
      call build_initial_graph(y,mval,affinity,status)
      if (status/=sgt_ok) then; result%status=status; return; end if
      n=size(affinity,1)
      if (kval<1 .or. kval>=n) then; result%status=sgt_invalid_input; return; end if
      allocate(s(n,n)); s=1.0_dp/real(n,dp)
      la=laplacian_from_directed_affinity(affinity)
      call symmetric_eigen_jacobi(la,eig,eigvec,status)
      f=eigvec(:,1:kval)
      allocate(lseq(maxit+1),times(maxit+1),pvec(n),row(n)); lseq=0.0_dp; times=0.0_dp
      lseq(1)=lmd; call cpu_time(start); nused=1; converged=.false.
      do iter=1,maxit
         v=pairwise_matrix_rownorm2(f)
         do i=1,n
            pvec=affinity(i,:)-0.5_dp*lmd*v(i,:)
            call project_simplex(pvec,row,status)
            s(i,:)=row
         end do
         ls=laplacian_from_directed_affinity(s)
         call symmetric_eigen_jacobi(ls,eig,eigvec,status)
         f=eigvec(:,1:kval)
         nzero=count(abs(eig)<eig_tol)
         nused=iter+1; call cpu_time(now); times(nused)=now-start
         if (nzero==kval) then
            converged=.true.; exit
         else if (nzero>kval) then
            lmd=0.5_dp*lmd
         else
            lmd=2.0_dp*lmd
         end if
         lseq(nused)=lmd
      end do
      where(abs(ls)<edge_tol) ls=0.0_dp
      result%laplacian=ls; result%adjacency=adjacency_from_laplacian(ls)
      result%eigenvalues=eig; result%eigenvectors=f
      result%parameter_history=lseq(1:nused); result%elapsed_time=times(1:nused)
      result%iterations=min(iter,maxit); result%convergence=converged
      if (converged) then; result%status=sgt_ok; else; result%status=sgt_no_convergence; end if
   end subroutine cluster_k_component_graph

   pure function degree_operator_matrix(p) result(s)
      integer, intent(in) :: p
      real(dp), allocatable :: s(:,:)
      integer :: m,i,j,k
      m=p*(p-1)/2; allocate(s(p,m)); s=0.0_dp; k=0
      do i=1,p-1
         do j=i+1,p
            k=k+1; s(i,k)=1.0_dp; s(j,k)=1.0_dp
         end do
      end do
   end function degree_operator_matrix

   subroutine learn_smooth_graph(x,result,alpha,beta,step_size,maxiter,tolerance)
      real(dp), intent(in) :: x(:,:)
      type(graph_result), intent(out) :: result
      real(dp), intent(in), optional :: alpha,beta,step_size,tolerance
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: sop(:,:),corr(:,:),corrinv(:,:),wk(:),dk(:),wkprev(:),dkprev(:)
      real(dp), allocatable :: yk(:),ykbar(:),pk(:),pkbar(:),qk(:),qkbar(:),z(:),dist(:,:)
      real(dp) :: alpha_v,beta_v,step_v,tol,mu,eps,gamma,denw,dend
      integer :: p,maxit,k,status
      logical :: converged
      result=graph_result()
      p=size(x,1)
      if (p<2 .or. size(x,2)<2) then; result%status=sgt_invalid_input; return; end if
      alpha_v=1e-2_dp; if (present(alpha)) alpha_v=alpha
      beta_v=1e-4_dp; if (present(beta)) beta_v=beta
      step_v=1e-2_dp; if (present(step_size)) step_v=step_size
      maxit=1000; if (present(maxiter)) maxit=max(1,maxiter)
      tol=1e-4_dp; if (present(tolerance)) tol=tolerance
      sop=degree_operator_matrix(p)
      call rows_correlation(x,corr,status)
      call symmetric_pseudoinverse(corr,corrinv,status)
      call initialize_weights_naive(corrinv,wk,status)
      dk=D(wk)
      mu=2.0_dp*beta_v+sqrt(2.0_dp*real(p-1,dp))
      eps=0.0_dp
      gamma=step_v*((1.0_dp-eps)/mu-eps)+eps
      dist=pairwise_matrix_rownorm2(x); z=upper_view_vec(dist)
      if (vector_norm2(z)>0.0_dp) z=z/vector_norm2(z)
      allocate(wkprev(size(wk)),dkprev(size(dk)),yk(size(wk)),ykbar(size(dk)), &
         pk(size(wk)),pkbar(size(dk)),qk(size(wk)),qkbar(size(dk)))
      converged=.false.
      do k=1,maxit
         wkprev=wk; dkprev=dk
         yk=wk-gamma*(2.0_dp*beta_v*wk+matmul(transpose(sop),dk))
         ykbar=dk+gamma*matmul(sop,wk)
         pk=max(0.0_dp,yk-2.0_dp*gamma*z)
         pkbar=0.5_dp*(ykbar-sqrt(ykbar*ykbar+4.0_dp*alpha_v*gamma))
         qk=pk-gamma*(2.0_dp*beta_v*pk+matmul(transpose(sop),pkbar))
         qkbar=pkbar+gamma*matmul(sop,pk)
         wk=max(wk-yk+qk,0.0_dp); dk=dk-ykbar+qkbar
         denw=max(vector_norm2(wkprev),tiny(1.0_dp)); dend=max(vector_norm2(dkprev),tiny(1.0_dp))
         if (vector_norm2(wk-wkprev)/denw<tol .and. vector_norm2(dk-dkprev)/dend<tol) then
            converged=.true.; exit
         end if
      end do
      result%weights=wk; result%laplacian=L(wk); result%adjacency=A(wk)
      result%iterations=min(k,maxit); result%convergence=converged
      if (converged) then; result%status=sgt_ok; else; result%status=sgt_no_convergence; end if
   end subroutine learn_smooth_graph

   subroutine update_signal_rep_l(y,alpha,beta,p,lw,status)
      real(dp), intent(in) :: y(:,:),alpha,beta
      integer, intent(in) :: p
      real(dp), allocatable, intent(out) :: lw(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: covy(:,:),linear(:),h(:,:),w(:),wnew(:),shifted(:)
      real(dp) :: floor_edge,total,lips,step,den
      integer :: m,k
      m=p*(p-1)/2
      covy=matmul(y,transpose(y))
      linear=alpha*Lstar(covy)
      h=2.0_dp*beta*Mmat(m)
      floor_edge=1e-6_dp/real(max(1,p-1),dp)
      total=0.5_dp*real(p,dp)-real(m,dp)*floor_edge
      if (total<=0.0_dp) then; allocate(lw(0,0)); status=sgt_invalid_input; return; end if
      allocate(w(m),wnew(m),shifted(m)); w=0.5_dp*real(p,dp)/real(m,dp)
      lips=maxval(sum(abs(h),dim=2)); step=1.0_dp/max(lips,tiny(1.0_dp))
      status=sgt_no_convergence
      do k=1,10000
         shifted=w-floor_edge-step*(matmul(h,w)+linear)
         call project_simplex_sum(shifted,total,wnew)
         wnew=wnew+floor_edge
         den=max(1.0_dp,vector_norm2(w))
         if (vector_norm2(wnew-w)/den<1e-10_dp) then; w=wnew; status=sgt_ok; exit; end if
         w=wnew
      end do
      lw=L(w)
   end subroutine update_signal_rep_l

   function signal_rep_objective(x,y,lw,alpha,beta) result(value)
      real(dp), intent(in) :: x(:,:),y(:,:),lw(:,:),alpha,beta
      real(dp) :: value
      real(dp) :: ly(size(y,1),size(y,2))
      ly=matmul(lw,y)
      value=frobenius_norm(x-y)**2+alpha*sum(y*ly)+beta*frobenius_norm(lw)**2
   end function signal_rep_objective

   subroutine learn_graph_sigrep(x,result,alpha,beta,maxiter,ftol)
      real(dp), intent(in) :: x(:,:)
      type(graph_result), intent(out) :: result
      real(dp), intent(in), optional :: alpha,beta,ftol
      integer, intent(in), optional :: maxiter
      real(dp), allocatable :: y(:,:),lw(:,:),system(:,:),ynew(:,:),obj(:)
      real(dp) :: alpha_v,beta_v,tol,fun0,fun,rel
      integer :: p,maxit,k,status,nused
      logical :: converged
      result=graph_result()
      p=size(x,1)
      if (p<2 .or. size(x,2)<1) then; result%status=sgt_invalid_input; return; end if
      alpha_v=1e-3_dp; if (present(alpha)) alpha_v=alpha
      beta_v=0.5_dp; if (present(beta)) beta_v=beta
      maxit=1000; if (present(maxiter)) maxit=max(1,maxiter)
      tol=1e-4_dp; if (present(ftol)) tol=ftol
      y=x; fun0=huge(1.0_dp); allocate(obj(maxit)); obj=0.0_dp; converged=.false.; nused=0
      do k=1,maxit
         call update_signal_rep_l(y,alpha_v,beta_v,p,lw,status)
         system=identity_matrix(p)+alpha_v*lw
         call solve_linear_system(system,x,ynew,status)
         if (status/=sgt_ok) exit
         y=ynew; fun=signal_rep_objective(x,y,lw,alpha_v,beta_v)
         nused=k; obj(k)=fun
         if (k>1) then
            rel=abs(fun-fun0)/max(abs(fun0),tiny(1.0_dp))
            if (rel<tol) then; converged=.true.; exit; end if
         end if
         fun0=fun
      end do
      result%laplacian=lw; result%adjacency=adjacency_from_laplacian(lw); result%smoothed_data=y
      if (nused>0) result%objective=obj(1:nused)
      result%iterations=min(k,maxit); result%convergence=converged
      if (converged) then; result%status=sgt_ok; else; result%status=sgt_no_convergence; end if
   end subroutine learn_graph_sigrep

   subroutine complete_adjacency_mask(p,mask)
      integer, intent(in) :: p
      real(dp), allocatable, intent(out) :: mask(:,:)
      integer :: i
      allocate(mask(p,p)); mask=1.0_dp
      do i=1,p; mask(i,i)=0.0_dp; end do
   end subroutine complete_adjacency_mask

   subroutine incidence_from_adjacency(mask,e)
      real(dp), intent(in) :: mask(:,:)
      real(dp), allocatable, intent(out) :: e(:,:)
      integer :: p,m,i,j,k
      p=size(mask,1); m=0
      do i=1,p-1; do j=i+1,p; if (mask(i,j)>0.0_dp) m=m+1; end do; end do
      allocate(e(p,m)); e=0.0_dp; k=0
      do i=1,p-1
         do j=i+1,p
            if (mask(i,j)>0.0_dp) then
               k=k+1; e(i,k)=1.0_dp; e(j,k)=-1.0_dp
            end if
         end do
      end do
   end subroutine incidence_from_adjacency

   subroutine learn_laplacian_gle_mm(s,result,a_mask,alpha,maxiter,reltol,record_objective)
      real(dp), intent(in) :: s(:,:)
      type(graph_result), intent(out) :: result
      real(dp), intent(in), optional :: a_mask(:,:),alpha,reltol
      integer, intent(in), optional :: maxiter
      logical, intent(in), optional :: record_objective
      real(dp), allocatable :: maskmat(:,:),sinv(:,:),allw(:),maskvec(:),w(:),wk(:),e(:,:),kmat(:,:),rmat(:,:)
      real(dp), allocatable :: g(:,:),gt(:,:),waug(:),gaug_t(:,:),gaug(:,:),middle(:,:),sol(:,:),q(:,:),z(:),obj(:)
      real(dp) :: alpha_v,tol,den
      integer :: p,m,maxit,k,status,i,nobj
      logical, allocatable :: active(:)
      logical :: rec,converged
      result=graph_result()
      p=size(s,1)
      if (p<2 .or. size(s,2)/=p) then; result%status=sgt_invalid_input; return; end if
      alpha_v=0.0_dp; if (present(alpha)) alpha_v=alpha
      maxit=10000; if (present(maxiter)) maxit=max(1,maxiter)
      tol=1e-5_dp; if (present(reltol)) tol=reltol
      rec=.false.; if (present(record_objective)) rec=record_objective
      if (present(a_mask)) then
         if (size(a_mask,1)/=p .or. size(a_mask,2)/=p) then; result%status=sgt_invalid_input; return; end if
         maskmat=a_mask
      else
         call complete_adjacency_mask(p,maskmat)
      end if
      maskvec=Ainv(maskmat); allocate(active(size(maskvec))); active=maskvec>0.0_dp
      call symmetric_pseudoinverse(s,sinv,status); call initialize_weights_naive(sinv,allw,status)
      w=pack(allw,active); wk=w; m=size(w)
      kmat=s+alpha_v*(identity_matrix(p)-1.0_dp)
      call incidence_from_adjacency(maskmat,e)
      rmat=matmul(transpose(e),matmul(kmat,e))
      allocate(g(p,m+1)); g(:,1:m)=e; g(:,m+1)=1.0_dp; gt=transpose(g)
      if (rec) then
         allocate(obj(maxit+1),z(size(allw)))
         obj=0.0_dp
         z=merge(allw,0.0_dp,active)
         obj(1)=vanilla_objective(L(z),kmat)
         nobj=1
      end if
      converged=.false.
      do k=1,maxit
         allocate(waug(m+1),gaug_t(m+1,p),gaug(p,m+1))
         waug(1:m)=wk; waug(m+1)=1.0_dp/real(p,dp)
         do i=1,m+1
            gaug_t(i,:)=gt(i,:)*waug(i)
         end do
         gaug=transpose(gaug_t)
         middle=matmul(gaug,transpose(g))
         call solve_linear_system(middle,gaug,sol,status)
         if (status/=sgt_ok) exit
         q=matmul(gaug_t,sol)
         do i=1,m
            if (rmat(i,i)>tiny(1.0_dp)) then; wk(i)=sqrt(max(q(i,i)/rmat(i,i),0.0_dp)); else; wk(i)=0.0_dp; end if
         end do
         if (rec) then
            z=0.0_dp; z=unpack(wk,active,0.0_dp); nobj=k+1; obj(nobj)=vanilla_objective(L(z),kmat)
         end if
         den=max(vector_norm2(w),tiny(1.0_dp)); converged=vector_norm2(w-wk)/den<tol
         deallocate(waug,gaug_t,gaug,middle,sol,q)
         if (converged .and. k>1) exit
         w=wk
      end do
      z=unpack(wk,active,0.0_dp)
      result%weights=z; result%laplacian=L(z); result%adjacency=A(z)
      result%iterations=min(k,maxit); result%convergence=converged
      if (rec) result%objective=obj(1:nobj)
      if (converged) then; result%status=sgt_ok; else; result%status=sgt_no_convergence; end if
   end subroutine learn_laplacian_gle_mm

   subroutine learn_laplacian_gle_admm(s,result,a_mask,alpha,rho,maxiter,reltol,record_objective)
      real(dp), intent(in) :: s(:,:)
      type(graph_result), intent(out) :: result
      real(dp), intent(in), optional :: a_mask(:,:),alpha,rho,reltol
      integer, intent(in), optional :: maxiter
      logical, intent(in), optional :: record_objective
      real(dp), allocatable :: maskmat(:,:),sinv(:,:),w(:),theta(:,:),thetak(:,:),yk(:,:),ck(:,:),cold(:,:)
      real(dp), allocatable :: kmat(:,:),pbase(:,:),gamma(:,:),eig(:),u(:,:),eval_map(:),xi(:,:),cktmp(:,:),rk(:,:),obj(:),times(:)
      real(dp) :: alpha_v,rho_v,tol,mu,tau,sres,rres,start,now,den
      integer :: p,maxit,k,status,nused,i,j
      logical :: rec,converged
      result=graph_result()
      p=size(s,1)
      if (p<2 .or. size(s,2)/=p) then; result%status=sgt_invalid_input; return; end if
      alpha_v=0.0_dp; if (present(alpha)) alpha_v=alpha
      rho_v=1.0_dp; if (present(rho)) rho_v=rho
      maxit=10000; if (present(maxiter)) maxit=max(1,maxiter)
      tol=1e-5_dp; if (present(reltol)) tol=reltol
      rec=.false.; if (present(record_objective)) rec=record_objective
      if (present(a_mask)) then; maskmat=a_mask; else; call complete_adjacency_mask(p,maskmat); end if
      call symmetric_pseudoinverse(s,sinv,status); call initialize_weights_naive(sinv,w,status)
      theta=L(w); yk=theta; ck=theta; cold=theta
      mu=2.0_dp; tau=2.0_dp; kmat=s+alpha_v*(identity_matrix(p)-1.0_dp)
      pbase=helmert_basis(p)
      allocate(times(maxit+1)); times=0.0_dp
      if (rec) then; allocate(obj(maxit+1)); obj=0.0_dp; obj(1)=vanilla_objective(theta,kmat); end if
      call cpu_time(start); nused=1; converged=.false.; thetak=theta
      do k=1,maxit
         gamma=matmul(transpose(pbase),matmul((kmat+yk)/rho_v-ck,pbase))
         call symmetric_eigen_jacobi(gamma,eig,u,status)
         allocate(eval_map(size(eig)))
         eval_map=0.5_dp*(sqrt(eig*eig+4.0_dp/rho_v)-eig)
         xi=spectral_reconstruction(u,eval_map); thetak=matmul(pbase,matmul(xi,transpose(pbase)))
         cktmp=yk/rho_v+thetak; ck=0.0_dp
         do i=1,p
            ck(i,i)=max(0.0_dp,cktmp(i,i))
            do j=i+1,p
               ck(i,j)=maskmat(i,j)*min(0.0_dp,cktmp(i,j)); ck(j,i)=ck(i,j)
            end do
         end do
         rk=thetak-ck; yk=yk+rho_v*rk
         nused=k+1; if (rec) obj(nused)=vanilla_objective(thetak,kmat)
         den=max(frobenius_norm(theta),tiny(1.0_dp)); converged=frobenius_norm(theta-thetak)/den<tol .and. k>1
         call cpu_time(now); times(nused)=now-start
         if (converged) exit
         sres=rho_v*frobenius_norm(cold-ck); rres=frobenius_norm(rk)
         if (rres>mu*sres) then; rho_v=rho_v*tau
         else if (sres>mu*rres) then; rho_v=rho_v/tau; end if
         theta=thetak; cold=ck; deallocate(eval_map)
      end do
      result%laplacian=thetak; result%adjacency=adjacency_from_laplacian(thetak)
      result%iterations=min(k,maxit); result%convergence=converged; result%elapsed_time=times(1:nused)
      if (rec) result%objective=obj(1:nused)
      if (converged) then; result%status=sgt_ok; else; result%status=sgt_no_convergence; end if
   end subroutine learn_laplacian_gle_admm

   subroutine update_sherman_morrison_diag(o,c,shift,idx)
      real(dp), intent(inout) :: o(:,:),c(:,:)
      real(dp), intent(in) :: shift
      integer, intent(in) :: idx
      real(dp) :: denom
      real(dp), allocatable :: col(:),row(:)
      o(idx,idx)=o(idx,idx)+shift
      denom=1.0_dp+shift*c(idx,idx)
      if (abs(denom)<=tiny(1.0_dp)) return
      col=c(:,idx); row=c(idx,:)
      c=c-(shift/denom)*outer_product(col,row)
   end subroutine update_sherman_morrison_diag

   subroutine learn_combinatorial_graph_laplacian(s,result,a_mask,alpha,reltol,max_cycle,regtype,record_objective)
      real(dp), intent(in) :: s(:,:)
      type(graph_result), intent(out) :: result
      real(dp), intent(in), optional :: a_mask(:,:),alpha,reltol
      integer, intent(in), optional :: max_cycle,regtype
      logical, intent(in), optional :: record_objective
      real(dp), allocatable :: maskmat(:,:),swork(:,:),hmat(:,:),kmat(:,:),o(:,:),c(:,:),oold(:,:),ou(:),cu(:)
      real(dp), allocatable :: ou_i(:,:),ku(:),b(:),beta_vec(:),hsub(:,:),bsub(:),x(:),ones(:),dshift(:),frob(:),times(:),obj(:)
      integer, allocatable :: minus(:),nzidx(:)
      real(dp) :: alpha_v,tol,dcvar,kuu,cuu_val,ouu,den,start,now
      integer :: n,maxcyc,reg,i,u,j,t,nnz,status,nused
      logical :: rec,shifting,converged
      result=graph_result()
      n=size(s,1)
      if (n<2 .or. size(s,2)/=n) then; result%status=sgt_invalid_input; return; end if
      alpha_v=0.0_dp; if (present(alpha)) alpha_v=alpha
      tol=1e-5_dp; if (present(reltol)) tol=reltol
      maxcyc=10000; if (present(max_cycle)) maxcyc=max(1,max_cycle)
      reg=1; if (present(regtype)) reg=regtype
      rec=.false.; if (present(record_objective)) rec=record_objective
      if (present(a_mask)) then; maskmat=a_mask; else; call complete_adjacency_mask(n,maskmat); end if
      swork=0.5_dp*(s+transpose(s)); allocate(ones(n)); ones=1.0_dp
      dcvar=dot_product(ones,matmul(swork,ones))/real(n*n,dp); shifting=abs(dcvar)<tol
      if (shifting) swork=swork+1.0_dp/real(n,dp)
      allocate(hmat(n,n))
      if (reg==1) then; hmat=alpha_v*(identity_matrix(n)-1.0_dp)
      else if (reg==2) then; hmat=alpha_v*(2.0_dp*identity_matrix(n)-1.0_dp)
      else; result%status=sgt_invalid_input; return; end if
      kmat=swork+hmat; allocate(o(n,n),c(n,n)); o=0.0_dp; c=0.0_dp
      do i=1,n
         if (kmat(i,i)<=tiny(1.0_dp)) then; result%status=sgt_invalid_input; return; end if
         o(i,i)=1.0_dp/kmat(i,i); c(i,i)=kmat(i,i)
      end do
      allocate(frob(maxcyc),times(maxcyc+1)); frob=0.0_dp; times=0.0_dp
      if (rec) then; allocate(obj(maxcyc+1)); obj=0.0_dp; obj(1)=vanilla_objective(o-1.0_dp/real(n,dp),kmat); end if
      call cpu_time(start); converged=.false.; nused=1
      do i=1,maxcyc
         oold=o
         do u=1,n
            allocate(minus(n-1)); minus=pack([(j,j=1,n)],[(j/=u,j=1,n)])
            kuu=kmat(u,u); allocate(ku(n-1),cu(n-1),ou_i(n-1,n-1),b(n-1),beta_vec(n-1))
            ku=kmat(minus,u); cu=c(minus,u); cuu_val=c(u,u)
            ou_i=c(minus,minus)-outer_product(cu,cu)/cuu_val
            b=ku/kuu+matmul(ou_i,ones(1:n-1))/real(n,dp); beta_vec=0.0_dp
            nnz=count(maskmat(minus,u)>0.0_dp)
            if (nnz>0) then
               allocate(nzidx(nnz)); nzidx=pack([(j,j=1,n-1)],maskmat(minus,u)>0.0_dp)
               allocate(hsub(nnz,nnz),bsub(nnz),x(nnz)); hsub=ou_i(nzidx,nzidx); bsub=b(nzidx); x=0.0_dp
               call nnqp_projected_gradient(hsub,bsub,x,status,max_iterations=5000,tolerance=1e-10_dp)
               beta_vec(nzidx)=-x
               deallocate(nzidx,hsub,bsub,x)
            end if
            ou=beta_vec+1.0_dp/real(n,dp)
            ouu=1.0_dp/kuu+dot_product(ou,matmul(ou_i,ou))
            o(u,u)=ouu; o(minus,u)=ou; o(u,minus)=ou
            den=ouu-dot_product(ou,matmul(ou_i,ou))
            cu=matmul(ou_i,ou)/den; cuu_val=1.0_dp/den
            c(u,u)=cuu_val; c(u,minus)=-cu; c(minus,u)=-cu
            c(minus,minus)=ou_i+outer_product(cu,cu)/cuu_val
            deallocate(minus,ku,cu,ou_i,b,beta_vec,ou)
         end do
         if (i>4) then
            dshift=matmul(o,ones)-1.0_dp
            do t=1,n
               if (abs(dshift(t))>1e-12_dp) call update_sherman_morrison_diag(o,c,-dshift(t),t)
            end do
         end if
         nused=i+1; if (rec) obj(nused)=vanilla_objective(o-1.0_dp/real(n,dp),kmat)
         call cpu_time(now); times(nused)=now-start
         frob(i)=frobenius_norm(oold-o)/max(frobenius_norm(oold),tiny(1.0_dp))
         if (i>6 .and. frob(i)<tol) then; converged=.true.; exit; end if
      end do
      result%laplacian=o-1.0_dp/real(n,dp); result%adjacency=adjacency_from_laplacian(result%laplacian)
      result%iterations=min(i,maxcyc); result%convergence=converged; result%elapsed_time=times(1:nused)
      result%parameter_history=frob(1:min(i,maxcyc)); if (rec) result%objective=obj(1:nused)
      if (converged) then; result%status=sgt_ok; else; result%status=sgt_no_convergence; end if
   end subroutine learn_combinatorial_graph_laplacian
end module sgt_other_learning
