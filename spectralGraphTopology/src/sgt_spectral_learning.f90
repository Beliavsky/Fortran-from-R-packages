! SPDX-License-Identifier: GPL-3.0-only
module sgt_spectral_learning
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use sgt_kinds, only : dp
   use sgt_status, only : sgt_ok, sgt_invalid_input, sgt_no_convergence
   use sgt_types, only : graph_result
   use sgt_linalg, only : symmetric_pseudoinverse, symmetric_eigen_jacobi, identity_matrix, &
      inverse_symmetric_positive_definite, is_positive_definite
   use sgt_operators, only : L, A
   use sgt_updates, only : initialize_weights_naive, initialize_weights_qp, &
      laplacian_w_update, joint_w_update, bipartite_w_update, laplacian_u_update, &
      bipartite_v_update, laplacian_lambda_update, bipartite_psi_update
   use sgt_objectives, only : laplacian_negative_log_likelihood, laplacian_prior, &
      bipartite_negative_log_likelihood, bipartite_prior, joint_prior, &
      bipartite_objective
   use sgt_initial_graph, only : build_initial_graph, laplacian_from_directed_affinity
   implicit none
   private
   public :: learn_k_component_graph, learn_cospectral_graph
   public :: learn_bipartite_graph, learn_bipartite_k_component_graph
contains
   subroutine prepare_input(s,is_data_matrix,m,swork,sinv,status)
      real(dp), intent(in) :: s(:,:)
      logical, intent(in) :: is_data_matrix
      integer, intent(in) :: m
      real(dp), allocatable, intent(out) :: swork(:,:),sinv(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: affinity(:,:),l0(:,:)
      logical :: data_mode
      data_mode=is_data_matrix .or. size(s,1)/=size(s,2)
      if (data_mode) then
         call build_initial_graph(s,m,affinity,status)
         if (status/=sgt_ok) then
            allocate(swork(0,0),sinv(0,0)); return
         end if
         l0=laplacian_from_directed_affinity(affinity)
         call symmetric_pseudoinverse(l0,swork,status)
         if (status/=sgt_ok .and. status/=sgt_no_convergence) then
            allocate(sinv(0,0)); return
         end if
         sinv=l0
      else
         if (size(s,1)<2) then
            allocate(swork(0,0),sinv(0,0)); status=sgt_invalid_input; return
         end if
         swork=0.5_dp*(s+transpose(s))
         call symmetric_pseudoinverse(swork,sinv,status)
      end if
   end subroutine prepare_input

   subroutine initial_weights(sinv,use_qp,w,status,floor_value)
      real(dp), intent(in) :: sinv(:,:)
      logical, intent(in) :: use_qp
      real(dp), allocatable, intent(out) :: w(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: floor_value
      if (use_qp) then
         call initialize_weights_qp(sinv,w,status)
         if (present(floor_value) .and. size(w)>0) w=w+max(0.0_dp,floor_value)
      else
         if (present(floor_value)) then
            call initialize_weights_naive(sinv,w,status,floor_value)
         else
            call initialize_weights_naive(sinv,w,status)
         end if
      end if
   end subroutine initial_weights

   subroutine learn_k_component_graph(s,result,is_data_matrix,k,use_qp,lb,ub,alpha,beta, &
      beta_max,fix_beta,rho,m,eps,maxiter,abstol,reltol,eigtol,record_objective,record_weights)
      real(dp), intent(in) :: s(:,:)
      type(graph_result), intent(out) :: result
      logical, intent(in), optional :: is_data_matrix,use_qp,fix_beta,record_objective,record_weights
      integer, intent(in), optional :: k,m,maxiter
      real(dp), intent(in), optional :: lb,ub,alpha,beta,beta_max,rho,eps,abstol,reltol,eigtol
      real(dp), allocatable :: sw(:,:),sinv(:,:),w0(:),w(:),lw0(:,:),lw(:,:),u0(:,:),u(:,:)
      real(dp), allocatable :: lambda0(:),lambda(:),hmat(:,:),kmat(:,:),eig(:),eigvec(:,:)
      real(dp), allocatable :: err(:),obj(:),ll(:),times(:),betaseq(:),whist(:,:)
      real(dp) :: lower,upper,alpha_v,beta_v,beta_limit,rho_v,eps_v,abs_v,rel_v,eig_v
      real(dp) :: start,now
      integer :: kval,mval,maxit,status,i,n,nzero,nused
      logical :: data_mode,qp_mode,fix_mode,rec_obj,rec_w,converged

      result=graph_result()
      data_mode=.false.; if (present(is_data_matrix)) data_mode=is_data_matrix
      qp_mode=.false.; if (present(use_qp)) qp_mode=use_qp
      fix_mode=.true.; if (present(fix_beta)) fix_mode=fix_beta
      rec_obj=.false.; if (present(record_objective)) rec_obj=record_objective
      rec_w=.false.; if (present(record_weights)) rec_w=record_weights
      allocate(whist(0,0))
      kval=1; if (present(k)) kval=k
      mval=7; if (present(m)) mval=m
      maxit=10000; if (present(maxiter)) maxit=max(1,maxiter)
      lower=0.0_dp; if (present(lb)) lower=lb
      upper=1e4_dp; if (present(ub)) upper=ub
      alpha_v=0.0_dp; if (present(alpha)) alpha_v=alpha
      beta_v=1e4_dp; if (present(beta)) beta_v=beta
      beta_limit=1e6_dp; if (present(beta_max)) beta_limit=beta_max
      rho_v=1e-2_dp; if (present(rho)) rho_v=rho
      eps_v=1e-4_dp; if (present(eps)) eps_v=eps
      abs_v=1e-6_dp; if (present(abstol)) abs_v=abstol
      rel_v=1e-4_dp; if (present(reltol)) rel_v=reltol
      eig_v=1e-9_dp; if (present(eigtol)) eig_v=eigtol
      if (beta_v<=0.0_dp .or. kval<1) then
         result%status=sgt_invalid_input; return
      end if
      call prepare_input(s,data_mode,mval,sw,sinv,status)
      if (status/=sgt_ok .and. status/=sgt_no_convergence) then
         result%status=status; return
      end if
      n=size(sw,1)
      if (kval>=n) then
         result%status=sgt_invalid_input; return
      end if
      call initial_weights(sinv,qp_mode,w0,status)
      lw0=L(w0)
      allocate(hmat(n,n)); hmat=alpha_v*(identity_matrix(n)-1.0_dp)
      kmat=sw+hmat
      call laplacian_u_update(lw0,kval,u0,status)
      call laplacian_lambda_update(lower,upper,beta_v,u0,lw0,kval,lambda0,status)
      allocate(times(maxit+1),betaseq(maxit+1)); times=0.0_dp; betaseq=0.0_dp
      betaseq(1)=beta_v
      if (rec_obj) then
         allocate(obj(maxit+1),ll(maxit+1)); obj=0.0_dp; ll=0.0_dp
         ll(1)=laplacian_negative_log_likelihood(lw0,lambda0,kmat)
         obj(1)=ll(1)+laplacian_prior(beta_v,lw0,lambda0,u0)
      end if
      if (rec_w) then
         deallocate(whist); allocate(whist(size(w0),maxit+1)); whist=0.0_dp; whist(:,1)=w0
      end if
      call cpu_time(start); converged=.false.; nused=1
      w=w0; lw=lw0; u=u0; lambda=lambda0
      do i=1,maxit
         call laplacian_w_update(w0,lw0,u0,beta_v,lambda0,kmat,w,status)
         lw=L(w)
         call laplacian_u_update(lw,kval,u,status)
         call laplacian_lambda_update(lower,upper,beta_v,u,lw,kval,lambda,status)
         nused=i+1
         if (rec_obj) then
            ll(nused)=laplacian_negative_log_likelihood(lw,lambda,kmat)
            obj(nused)=ll(nused)+laplacian_prior(beta_v,lw,lambda,u)
         end if
         if (rec_w) whist(:,nused)=w
         err=abs(w0-w)
         converged=all(err<=0.5_dp*rel_v*(w+w0)) .or. all(err<=abs_v)
         call cpu_time(now); times(nused)=now-start
         if (.not.fix_mode) then
            call symmetric_eigen_jacobi(lw,eig,eigvec,status)
            nzero=count(abs(eig)<eig_v)
            if (kval<=nzero) then
               beta_v=(1.0_dp+rho_v)*beta_v
            else
               beta_v=beta_v/(1.0_dp+rho_v)
            end if
            beta_v=min(beta_v,beta_limit)
         end if
         betaseq(nused)=beta_v
         if (converged) exit
         w0=w; lw0=lw; u0=u; lambda0=lambda
         kmat=sw+hmat/(-lw+eps_v)
      end do
      result%laplacian=lw; result%adjacency=A(w); result%weights=w
      result%eigenvalues=lambda; result%eigenvectors=u
      result%iterations=min(i,maxit); result%convergence=converged; result%beta=beta_v
      result%elapsed_time=times(1:nused); result%parameter_history=betaseq(1:nused)
      if (rec_obj) then
         result%objective=obj(1:nused); result%negative_log_likelihood=ll(1:nused)
      end if
      if (rec_w) result%weight_history=whist(:,1:nused)
      if (converged) then
         result%status=sgt_ok
      else
         result%status=sgt_no_convergence
      end if
   end subroutine learn_k_component_graph

   subroutine learn_cospectral_graph(s,fixed_lambda,result,is_data_matrix,k,use_qp,alpha,beta, &
      beta_max,fix_beta,rho,m,maxiter,abstol,reltol,eigtol,record_objective,record_weights)
      real(dp), intent(in) :: s(:,:),fixed_lambda(:)
      type(graph_result), intent(out) :: result
      logical, intent(in), optional :: is_data_matrix,use_qp,fix_beta,record_objective,record_weights
      integer, intent(in), optional :: k,m,maxiter
      real(dp), intent(in), optional :: alpha,beta,beta_max,rho,abstol,reltol,eigtol
      real(dp), allocatable :: sw(:,:),sinv(:,:),w0(:),w(:),lw0(:,:),lw(:,:),u0(:,:),u(:,:)
      real(dp), allocatable :: hmat(:,:),kmat(:,:),eig(:),eigvec(:,:),err(:),obj(:),ll(:),times(:),betaseq(:),whist(:,:)
      real(dp) :: alpha_v,beta_v,beta_limit,rho_v,abs_v,rel_v,eig_v,start,now
      integer :: kval,mval,maxit,status,i,n,nzero,nused
      logical :: data_mode,qp_mode,fix_mode,rec_obj,rec_w,converged

      result=graph_result()
      data_mode=.false.; if (present(is_data_matrix)) data_mode=is_data_matrix
      qp_mode=.false.; if (present(use_qp)) qp_mode=use_qp
      fix_mode=.true.; if (present(fix_beta)) fix_mode=fix_beta
      rec_obj=.false.; if (present(record_objective)) rec_obj=record_objective
      rec_w=.false.; if (present(record_weights)) rec_w=record_weights
      allocate(whist(0,0))
      kval=1; if (present(k)) kval=k
      mval=7; if (present(m)) mval=m
      maxit=10000; if (present(maxiter)) maxit=max(1,maxiter)
      alpha_v=0.0_dp; if (present(alpha)) alpha_v=alpha
      beta_v=1e4_dp; if (present(beta)) beta_v=beta
      beta_limit=1e6_dp; if (present(beta_max)) beta_limit=beta_max
      rho_v=1e-2_dp; if (present(rho)) rho_v=rho
      abs_v=1e-6_dp; if (present(abstol)) abs_v=abstol
      rel_v=1e-4_dp; if (present(reltol)) rel_v=reltol
      eig_v=1e-9_dp; if (present(eigtol)) eig_v=eigtol
      call prepare_input(s,data_mode,mval,sw,sinv,status)
      if (status/=sgt_ok .and. status/=sgt_no_convergence) then
         result%status=status; return
      end if
      n=size(sw,1)
      if (kval<1 .or. kval>=n .or. size(fixed_lambda)/=n-kval .or. minval(fixed_lambda)<=0.0_dp) then
         result%status=sgt_invalid_input; return
      end if
      call initial_weights(sinv,qp_mode,w0,status); lw0=L(w0)
      allocate(hmat(n,n)); hmat=alpha_v*(2.0_dp*identity_matrix(n)-1.0_dp); kmat=sw+hmat
      call laplacian_u_update(lw0,kval,u0,status)
      allocate(times(maxit+1),betaseq(maxit+1)); times=0.0_dp; betaseq=0.0_dp; betaseq(1)=beta_v
      if (rec_obj) then
         allocate(obj(maxit+1),ll(maxit+1)); obj=0.0_dp; ll=0.0_dp
         ll(1)=laplacian_negative_log_likelihood(lw0,fixed_lambda,kmat)
         obj(1)=ll(1)+laplacian_prior(beta_v,lw0,fixed_lambda,u0)
      end if
      if (rec_w) then
         deallocate(whist); allocate(whist(size(w0),maxit+1)); whist=0.0_dp; whist(:,1)=w0
      end if
      call cpu_time(start); converged=.false.; nused=1
      w=w0; lw=lw0; u=u0
      do i=1,maxit
         call laplacian_w_update(w0,lw0,u0,beta_v,fixed_lambda,kmat,w,status)
         lw=L(w); call laplacian_u_update(lw,kval,u,status)
         nused=i+1
         if (rec_obj) then
            ll(nused)=laplacian_negative_log_likelihood(lw,fixed_lambda,kmat)
            obj(nused)=ll(nused)+laplacian_prior(beta_v,lw,fixed_lambda,u)
         end if
         if (rec_w) whist(:,nused)=w
         err=abs(w0-w); converged=all(err<=0.5_dp*rel_v*(w+w0)) .or. all(err<=abs_v)
         call cpu_time(now); times(nused)=now-start
         if (.not.fix_mode) then
            call symmetric_eigen_jacobi(lw,eig,eigvec,status); nzero=count(abs(eig)<eig_v)
            if (kval<=nzero) then; beta_v=(1.0_dp+rho_v)*beta_v
            else; beta_v=beta_v/(1.0_dp+rho_v); end if
            beta_v=min(beta_v,beta_limit)
         end if
         betaseq(nused)=beta_v
         if (converged) exit
         w0=w; lw0=lw; u0=u
      end do
      result%laplacian=lw; result%adjacency=A(w); result%weights=w
      result%eigenvalues=fixed_lambda; result%eigenvectors=u; result%iterations=min(i,maxit)
      result%convergence=converged; result%beta=beta_v; result%elapsed_time=times(1:nused)
      result%parameter_history=betaseq(1:nused)
      if (rec_obj) then; result%objective=obj(1:nused); result%negative_log_likelihood=ll(1:nused); end if
      if (rec_w) result%weight_history=whist(:,1:nused)
      if (converged) then; result%status=sgt_ok; else; result%status=sgt_no_convergence; end if
   end subroutine learn_cospectral_graph

   subroutine learn_bipartite_graph(s,result,is_data_matrix,z,nu,alpha,use_qp,m,maxiter, &
      abstol,reltol,record_weights,record_objective)
      real(dp), intent(in) :: s(:,:)
      type(graph_result), intent(out) :: result
      logical, intent(in), optional :: is_data_matrix,use_qp,record_weights,record_objective
      integer, intent(in), optional :: z,m,maxiter
      real(dp), intent(in), optional :: nu,alpha,abstol,reltol
      real(dp), allocatable :: sw(:,:),sinv(:,:),w0(:),w(:),aw0(:,:),aw(:,:),lw(:,:),v0(:,:),v(:,:),psi0(:),psi(:)
      real(dp), allocatable :: jmat(:,:),hmat(:,:),kmat(:,:),eig(:),eigvec(:,:),err(:),obj(:),ll(:),times(:),lipseq(:),whist(:,:)
      real(dp) :: nu_v,alpha_v,abs_v,rel_v,lips,fun0,funt,start,now
      integer :: zval,mval,maxit,status,i,n,nused,backtrack
      logical :: data_mode,qp_mode,rec_w,rec_obj,converged,accepted

      result=graph_result()
      data_mode=.false.; if (present(is_data_matrix)) data_mode=is_data_matrix
      qp_mode=.false.; if (present(use_qp)) qp_mode=use_qp
      rec_w=.false.; if (present(record_weights)) rec_w=record_weights
      allocate(whist(0,0))
      rec_obj=.true.; if (present(record_objective)) rec_obj=record_objective
      zval=0; if (present(z)) zval=z
      mval=7; if (present(m)) mval=m
      maxit=10000; if (present(maxiter)) maxit=max(1,maxiter)
      nu_v=1e4_dp; if (present(nu)) nu_v=nu
      alpha_v=0.0_dp; if (present(alpha)) alpha_v=alpha
      abs_v=1e-6_dp; if (present(abstol)) abs_v=abstol
      rel_v=1e-4_dp; if (present(reltol)) rel_v=reltol
      if (nu_v<=0.0_dp) then; result%status=sgt_invalid_input; return; end if
      call prepare_input(s,data_mode,mval,sw,sinv,status)
      if (status/=sgt_ok .and. status/=sgt_no_convergence) then; result%status=status; return; end if
      n=size(sw,1)
      if (zval<0 .or. zval>=n .or. mod(n-zval,2)/=0) then; result%status=sgt_invalid_input; return; end if
      allocate(jmat(n,n),hmat(n,n)); jmat=1.0_dp/real(n,dp)
      hmat=alpha_v*(2.0_dp*identity_matrix(n)-1.0_dp); kmat=sw+hmat
      call initial_weights(sinv,qp_mode,w0,status,1e-4_dp)
      call symmetric_eigen_jacobi(L(w0)+jmat,eig,eigvec,status)
      lips=1.0_dp/max(minval(eig),1e-12_dp)
      aw0=A(w0); call bipartite_v_update(aw0,zval,v0,status); call bipartite_psi_update(v0,aw0,psi0,status)
      allocate(times(maxit+1),lipseq(20*maxit+1)); times=0.0_dp; lipseq=0.0_dp; lipseq(1)=lips
      allocate(obj(maxit+1),ll(maxit+1)); obj=0.0_dp; ll=0.0_dp
      ll(1)=bipartite_negative_log_likelihood(L(w0),kmat,jmat)
      obj(1)=ll(1)+bipartite_prior(nu_v,aw0,psi0,v0); fun0=obj(1)
      if (rec_w) then; deallocate(whist); allocate(whist(size(w0),maxit+1)); whist=0.0_dp; whist(:,1)=w0; end if
      call cpu_time(start); nused=1; backtrack=1; converged=.false.; w=w0; aw=aw0; v=v0; psi=psi0
      do i=1,maxit
         accepted=.false.
         do while (.not.accepted)
            call bipartite_w_update(w0,aw0,v0,nu_v,psi0,kmat,jmat,lips,w,status)
            lw=L(w); aw=A(w)
            if (status==sgt_ok) then
               if (is_positive_definite(lw+jmat)) then
                  funt=bipartite_objective(aw,lw,v0,psi0,kmat,jmat,nu_v)
               else
                  funt=huge(1.0_dp)
               end if
            else
               funt=huge(1.0_dp)
            end if
            if (backtrack<=size(lipseq)) then; lipseq(backtrack)=lips; backtrack=backtrack+1; end if
            if (.not.ieee_is_finite(funt) .or. fun0<funt) then
               lips=2.0_dp*lips
               if (lips>1e20_dp) exit
            else
               lips=max(0.5_dp*lips,1e-12_dp); accepted=.true.
            end if
         end do
         if (.not.accepted) exit
         call bipartite_v_update(aw,zval,v,status); call bipartite_psi_update(v,aw,psi,status)
         nused=i+1
         ll(nused)=bipartite_negative_log_likelihood(lw,kmat,jmat)
         obj(nused)=ll(nused)+bipartite_prior(nu_v,aw,psi,v)
         if (rec_w) whist(:,nused)=w
         err=abs(w0-w); converged=all(err<=0.5_dp*rel_v*(w+w0)) .or. all(err<=abs_v)
         call cpu_time(now); times(nused)=now-start
         if (converged) exit
         fun0=obj(nused); w0=w; aw0=aw; v0=v; psi0=psi
      end do
      result%laplacian=L(w); result%adjacency=A(w); result%weights=w
      result%auxiliary_eigenvalues=psi; result%auxiliary_eigenvectors=v
      result%iterations=min(i,maxit); result%convergence=converged; result%nu=nu_v; result%lipschitz=lips
      result%elapsed_time=times(1:nused); result%parameter_history=lipseq(1:max(1,backtrack-1))
      if (rec_obj) then; result%objective=obj(1:nused); result%negative_log_likelihood=ll(1:nused); end if
      if (rec_w) result%weight_history=whist(:,1:nused)
      if (converged) then; result%status=sgt_ok; else; result%status=sgt_no_convergence; end if
   end subroutine learn_bipartite_graph

   subroutine learn_bipartite_k_component_graph(s,result,is_data_matrix,z,k,use_qp,m,alpha,beta, &
      rho,fix_beta,beta_max,nu,lb,ub,maxiter,abstol,reltol,eigtol,record_weights,record_objective)
      real(dp), intent(in) :: s(:,:)
      type(graph_result), intent(out) :: result
      logical, intent(in), optional :: is_data_matrix,use_qp,fix_beta,record_weights,record_objective
      integer, intent(in), optional :: z,k,m,maxiter
      real(dp), intent(in), optional :: alpha,beta,rho,beta_max,nu,lb,ub,abstol,reltol,eigtol
      real(dp), allocatable :: sw(:,:),sinv(:,:),w0(:),w(:),lw0(:,:),lw(:,:),aw0(:,:),aw(:,:)
      real(dp), allocatable :: u0(:,:),u(:,:),v0(:,:),v(:,:),lambda0(:),lambda(:),psi0(:),psi(:)
      real(dp), allocatable :: hmat(:,:),kmat(:,:),eig(:),eigvec(:,:),err(:),obj(:),ll(:),times(:),betaseq(:),whist(:,:)
      real(dp) :: alpha_v,beta_v,rho_v,beta_limit,nu_v,lower,upper,abs_v,rel_v,eig_v,start,now
      integer :: zval,kval,mval,maxit,status,i,n,nused,nzero
      logical :: data_mode,qp_mode,fix_mode,rec_w,rec_obj,converged

      result=graph_result()
      data_mode=.false.; if (present(is_data_matrix)) data_mode=is_data_matrix
      qp_mode=.false.; if (present(use_qp)) qp_mode=use_qp
      fix_mode=.true.; if (present(fix_beta)) fix_mode=fix_beta
      rec_w=.false.; if (present(record_weights)) rec_w=record_weights
      allocate(whist(0,0))
      rec_obj=.false.; if (present(record_objective)) rec_obj=record_objective
      zval=0; if (present(z)) zval=z
      kval=1; if (present(k)) kval=k
      mval=7; if (present(m)) mval=m
      maxit=10000; if (present(maxiter)) maxit=max(1,maxiter)
      alpha_v=0.0_dp; if (present(alpha)) alpha_v=alpha
      beta_v=1e4_dp; if (present(beta)) beta_v=beta
      rho_v=1e-2_dp; if (present(rho)) rho_v=rho
      beta_limit=1e6_dp; if (present(beta_max)) beta_limit=beta_max
      nu_v=1e4_dp; if (present(nu)) nu_v=nu
      lower=0.0_dp; if (present(lb)) lower=lb
      upper=1e4_dp; if (present(ub)) upper=ub
      abs_v=1e-6_dp; if (present(abstol)) abs_v=abstol
      rel_v=1e-4_dp; if (present(reltol)) rel_v=reltol
      eig_v=1e-9_dp; if (present(eigtol)) eig_v=eigtol
      if (beta_v<0.0_dp .or. nu_v<0.0_dp .or. beta_v+nu_v<=0.0_dp) then; result%status=sgt_invalid_input; return; end if
      call prepare_input(s,data_mode,mval,sw,sinv,status)
      if (status/=sgt_ok .and. status/=sgt_no_convergence) then; result%status=status; return; end if
      n=size(sw,1)
      if (kval<1 .or. kval>=n .or. zval<0 .or. zval>=n .or. mod(n-zval,2)/=0) then
         result%status=sgt_invalid_input; return
      end if
      allocate(hmat(n,n)); hmat=alpha_v*(2.0_dp*identity_matrix(n)-1.0_dp); kmat=sw+hmat
      call initial_weights(sinv,qp_mode,w0,status)
      lw0=L(w0); aw0=A(w0)
      call bipartite_v_update(aw0,zval,v0,status); call bipartite_psi_update(v0,aw0,psi0,status)
      call laplacian_u_update(lw0,kval,u0,status)
      call laplacian_lambda_update(lower,upper,max(beta_v,epsilon(1.0_dp)),u0,lw0,kval,lambda0,status)
      allocate(times(maxit+1),betaseq(maxit+1)); times=0.0_dp; betaseq=0.0_dp; betaseq(1)=beta_v
      if (rec_obj) then
         allocate(obj(maxit+1),ll(maxit+1)); obj=0.0_dp; ll=0.0_dp
         ll(1)=laplacian_negative_log_likelihood(lw0,lambda0,kmat)
         obj(1)=ll(1)+joint_prior(beta_v,nu_v,lw0,aw0,u0,v0,lambda0,psi0)
      end if
      if (rec_w) then; deallocate(whist); allocate(whist(size(w0),maxit+1)); whist=0.0_dp; whist(:,1)=w0; end if
      call cpu_time(start); converged=.false.; nused=1
      w=w0; lw=lw0; aw=aw0; u=u0; v=v0; lambda=lambda0; psi=psi0
      do i=1,maxit
         call joint_w_update(w0,lw0,aw0,u0,v0,lambda0,psi0,beta_v,nu_v,kmat,w,status)
         lw=L(w); aw=A(w)
         call laplacian_u_update(lw,kval,u,status); call bipartite_v_update(aw,zval,v,status)
         if (beta_v>0.0_dp) then
            call laplacian_lambda_update(lower,upper,beta_v,u,lw,kval,lambda,status)
         else
            lambda=lambda0
         end if
         call bipartite_psi_update(v,aw,psi,status)
         nused=i+1
         if (rec_obj) then
            ll(nused)=laplacian_negative_log_likelihood(lw,lambda,kmat)
            obj(nused)=ll(nused)+joint_prior(beta_v,nu_v,lw,aw,u,v,lambda,psi)
         end if
         if (rec_w) whist(:,nused)=w
         err=abs(w0-w); converged=all(err<=0.5_dp*rel_v*(w+w0)) .or. all(err<=abs_v)
         call cpu_time(now); times(nused)=now-start
         call symmetric_eigen_jacobi(lw,eig,eigvec,status)
         if (.not.fix_mode) then
            nzero=count(abs(eig)<eig_v)
            if (kval<nzero) then; beta_v=(1.0_dp+rho_v)*beta_v
            else if (kval>nzero) then; beta_v=beta_v/(1.0_dp+rho_v); end if
            beta_v=min(beta_v,beta_limit)
         end if
         betaseq(nused)=beta_v
         if (converged) exit
         w0=w; lw0=lw; aw0=aw; u0=u; v0=v; lambda0=lambda; psi0=psi
      end do
      result%laplacian=lw; result%adjacency=aw; result%weights=w
      result%eigenvalues=lambda; result%eigenvectors=u
      result%auxiliary_eigenvalues=psi; result%auxiliary_eigenvectors=v
      result%iterations=min(i,maxit); result%convergence=converged; result%beta=beta_v; result%nu=nu_v
      result%elapsed_time=times(1:nused); result%parameter_history=betaseq(1:nused)
      if (rec_obj) then; result%objective=obj(1:nused); result%negative_log_likelihood=ll(1:nused); end if
      if (rec_w) result%weight_history=whist(:,1:nused)
      if (converged) then; result%status=sgt_ok; else; result%status=sgt_no_convergence; end if
   end subroutine learn_bipartite_k_component_graph
end module sgt_spectral_learning
